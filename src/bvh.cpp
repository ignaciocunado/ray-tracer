#include "bvh.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "extra.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>


// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
        glm::vec3 calcP = glm::vec3 { hitInfo.barycentricCoord.x * v0.position.x, hitInfo.barycentricCoord.x * v0.position.y, hitInfo.barycentricCoord.x * v0.position.z} + 
        glm::vec3 { hitInfo.barycentricCoord.y * v1.position.x, hitInfo.barycentricCoord.y * v1.position.y, hitInfo.barycentricCoord.y * v1.position.z } +
        glm::vec3 { hitInfo.barycentricCoord.z * v2.position.x, hitInfo.barycentricCoord.z * v2.position.y, hitInfo.barycentricCoord.z * v2.position.z };
        drawLine(calcP, hitInfo.normal, glm::vec3{0,0,1});
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex);

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    glm::vec3 minCoordinates, maxCoordinates; // minCoordinate[0] is minimum x-coordinate, 
                                              // minCoordinate[1] is minimum y-coordinate, and so on.

    // loop through each coordinate (x, then y, then z)
    for (int coord = 0; coord < 3; coord++) {
        float v0Coord = primitive.v0.position[coord];
        float v1Coord = primitive.v1.position[coord];
        float v2Coord = primitive.v2.position[coord];

        minCoordinates[coord] = std::min(v0Coord, std::min(v1Coord, v2Coord));
        maxCoordinates[coord] = std::max(v0Coord, std::max(v1Coord, v2Coord));
    }

    return { .lower = minCoordinates, .upper = maxCoordinates };
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    if (primitives.size() == 0) {
        return AxisAlignedBox {};
    }

    const AxisAlignedBox box0 = computePrimitiveAABB(primitives[0]);
    glm::vec3 minCoordinates = box0.lower;
    glm::vec3 maxCoordinates = box0.upper;

    for (int i = 1; i < primitives.size(); i++) {
        // Updating minCoordinates and maxCoordinates, so that primitives[i] is inside them.

        const AxisAlignedBox box = computePrimitiveAABB(primitives[i]);

        for (int coord = 0; coord < 3; coord++) {
            minCoordinates[coord] = std::min(minCoordinates[coord], box.lower[coord]);
            maxCoordinates[coord] = std::max(maxCoordinates[coord], box.upper[coord]);
        }
    }

    return { .lower = minCoordinates, .upper = maxCoordinates };
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    // Definition of the Centroid (of any figure) can be found at the top of this page:
    // https://en.wikipedia.org/wiki/Centroid

    return (primitive.v0.position + primitive.v1.position + primitive.v2.position) / 3.0f;
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    float xLength = aabb.upper.x - aabb.lower.x;
    float yLength = aabb.upper.y - aabb.lower.y;
    float zLength = aabb.upper.z - aabb.lower.z;

    if (xLength >= yLength && xLength >= zLength) {
        return 0;
    } else if (yLength >= zLength) {
        return 1;
    } else {
        return 2;
    }
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element 
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    using Primitive = BVHInterface::Primitive;

    // Take centroids of primitives and sort them by given axis.
    // Sorting is done with lambda function to avoid writing custom comparators.
    // Examples can be found in C++ reference: https://en.cppreference.com/w/cpp/algorithm/ranges/sort
    std::ranges::sort(primitives, [axis](const Primitive& a, const Primitive& b) {
        const glm::vec3 aCentroid = computePrimitiveCentroid(a);
        const glm::vec3 bCentroid = computePrimitiveCentroid(b);
        return aCentroid[axis] < bCentroid[axis];
    });

    // Calculate split position and return it.
    return (primitives.size() + 1) / 2;
}

// Helper function for ray-primitive (ray-triangle) intersection.
// Also updates hitInfo.
// Called by intersectRayWithBVH(...) and by intersectRayWithBVHWhenEnabledAccel(...).
bool intersectRayWithPrimitive(RenderState& state, const BVHInterface::Primitive& primitive, Ray& ray, HitInfo& hitInfo) {
    if (intersectRayWithTriangle(
            primitive.v0.position,
            primitive.v1.position,
            primitive.v2.position,
            ray, hitInfo)) {
        updateHitInfo(state, primitive, ray, hitInfo);
        return true;
    }

    return false;
}

// Helper function for traversing the BVH.
// Called by intersectRayWithBVH(...).
// Note that when calling this function it is assumed
// that state.features.enableAccelStructure is true.
bool intersectRayWithBVHWhenEnabledAccel(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Firstly, let's check if the ray would hit the AABB of the root.
    // Copy current ray. Make its t infinite (in case the ray.origin was inside AABB).
    Ray rayThisBox = ray;
    rayThisBox.t = std::numeric_limits<float>::max();
    if (!intersectRayWithShape(nodes[BVH::RootIndex].aabb, rayThisBox)) {
        // The ray does not intersect with the AABB of the root.
        // (Root's AABB contains ALL the primitives.)
        return false;
    }

    // Return value.
    bool isHit = false;

    // Indices of nodes to traverse.
    std::vector<uint32_t> nodeStack { BVH::RootIndex }; 

    while (nodeStack.size() > 0) {
        // Get the current node index from top of the stack, and pop it.
        int nodeIndex = nodeStack.back();
        nodeStack.pop_back();

        if (nodes[nodeIndex].isLeaf()) {
            for (uint32_t primitiveIndex = nodes[nodeIndex].primitiveOffset();
                 primitiveIndex < nodes[nodeIndex].primitiveOffset() + nodes[nodeIndex].primitiveCount();
                 primitiveIndex++) {
                isHit |= intersectRayWithPrimitive(state, primitives[primitiveIndex], ray, hitInfo);
            }
            continue;
        }

        // Current node is NOT a leaf, so has two children.
        
        // Copy current ray. Used to check intersection with children node AABBs.
        Ray rayLeftBox = ray;
        Ray rayRightBox = ray;
            
        // We want to set ray.t back to infinity since there might be the case when ray.origin is inside AABB.
        rayLeftBox.t = std::numeric_limits<float>::max();
        rayRightBox.t = std::numeric_limits<float>::max();

        uint32_t leftChildIndex = nodes[nodeIndex].leftChild();
        uint32_t rightChildIndex = nodes[nodeIndex].rightChild();

        bool hitLeftChild = intersectRayWithShape(nodes[leftChildIndex].aabb, rayLeftBox);
        bool hitRightChild = intersectRayWithShape(nodes[rightChildIndex].aabb, rayRightBox);

        if (hitLeftChild && hitRightChild) {
            // Push both children to nodeStack.
            // Firstly push the one that has ray.t larger, so that the one with
            // lower ray.t is on top of the stack and is traversed first.
            if (rayLeftBox.t < rayRightBox.t) {
                nodeStack.push_back(rightChildIndex);
                nodeStack.push_back(leftChildIndex);
            } else {
                nodeStack.push_back(leftChildIndex);
                nodeStack.push_back(rightChildIndex);
            }
        } else if (hitLeftChild) {
            nodeStack.push_back(leftChildIndex);
        } else if (hitRightChild) {
            nodeStack.push_back(rightChildIndex);
        }
    }
    
    return isHit;
}

// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default implementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Return value
    bool isHit = false;

    if (state.features.enableAccelStructure) {
        isHit |= intersectRayWithBVHWhenEnabledAccel(state, bvh, ray, hitInfo);
    } else {
        // Acceleration structures are NOT enabled.
        // Naive implementation; simply iterates over all primitives.
        std::span<const BVHInterface::Primitive> primitives = bvh.primitives();
        for (const auto& prim : primitives) {
            isHit |= intersectRayWithPrimitive(state, prim, ray, hitInfo);
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        isHit |= intersectRayWithShape(sphere, ray, hitInfo);

    return isHit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    uint32_t offset = static_cast<uint32_t>(m_primitives.size());
    uint32_t count = static_cast<uint32_t>(primitives.size());

    Node node {
        .aabb = aabb,
        .data = {
            // according to bvh_interface.h
            offset | Node::LeafBit,
            count
        }
    };

    // Copy the current set of primitives to the back of the primitives vector
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));

    return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    return Node {
        .aabb = aabb,
        .data = {
            // according to bvh_interface.h
            leftChildIndex,
            rightChildIndex }
    };
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    AxisAlignedBox aabb = computeSpanAABB(primitives);

    if (primitives.size() <= LeafSize) {
        // Current node is a leaf
        m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
        return;
    }

    // Current node is a node (NOT a leaf)

    size_t splitPosition = 0;

    if (features.extra.enableBvhSahBinning) {
        // TODO: this will be implemented in Extra features.
        // here, the primitives should be split into two subranges based on sah-binning
    } else {
        // Sort and split by median along the longest axis.
        uint32_t longestAxis = computeAABBLongestAxis(aabb);
        splitPosition = splitPrimitivesByMedian(aabb, longestAxis, primitives);
    }

    // Allocate left and right child nodes.
    uint32_t leftChildIndex = nextNodeIdx();
    uint32_t rightChildIndex = nextNodeIdx();

    // Fill in data for the current node.
    m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftChildIndex, rightChildIndex);

    // Build left and right children recursively.
    std::span<Primitive> leftChildPrimitives = primitives.subspan(0, splitPosition);
    std::span<Primitive> rightChildPrimitives = primitives.subspan(splitPosition, primitives.size() - splitPosition);
    buildRecursive(scene, features, leftChildPrimitives, leftChildIndex);
    buildRecursive(scene, features, rightChildPrimitives, rightChildIndex);
}

// Helper function which calculates the level (depth) for each node.
// levelsToFillIn - this vector will be filled with level for each node index.
// levelsToFillIn should be already of size of m_nodes.size() when passed to this function.
void BVH::calculateLevels(std::vector<uint32_t>& levelsToFillIn) {
    levelsToFillIn[RootIndex] = 0;

    for (int nodeIndex = 0; nodeIndex < m_nodes.size(); nodeIndex++) {
        if (nodeIndex == 1) {
            // This is a dummy node created in BVH constructor
            continue;
        }

        m_numLevels = std::max(m_numLevels, levelsToFillIn[nodeIndex]);
        const BVHInterface::Node& node = m_nodes[nodeIndex];
        if (!node.isLeaf()) {
            // Calculate depths of left and right children of the current node.
            levelsToFillIn[node.leftChild()] = levelsToFillIn[nodeIndex] + 1;
            levelsToFillIn[node.rightChild()] = levelsToFillIn[nodeIndex] + 1;
        }
    }
}

// TODO: Standard feature, or part of it
// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
    // Storing levels of each of the nodes.
    std::vector<uint32_t> levels(m_nodes.size());
    calculateLevels(levels);

    m_numLevels = 1;

    // For loop starts from 2. Node with index 1 is skipped since it was
    // added as a "dummy node" in the BVH constructor.
    for (int nodeIndex = 2; nodeIndex < m_nodes.size(); nodeIndex++) {
        m_numLevels = std::max(m_numLevels, levels[nodeIndex] + 1);
    }
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    m_numLeaves = 0;
    for (int i = 0; i < m_nodes.size(); i++) {
        // Check if current node is leaf, and is not a dummy node (index 1),
        // which was added in BVH constructor.
        if (i != 1 && m_nodes[i].isLeaf()) {
            m_numLeaves++;
        }
    }
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
    // Storing levels of each of the nodes.
    std::vector<uint32_t> levels(m_nodes.size());
    calculateLevels(levels);

    const glm::vec3 color = glm::vec3(0.5f, 1.0f, 0.5f);
    const float transparency = 0.6f;

    for (int nodeIndex = 0; nodeIndex < m_nodes.size(); nodeIndex++) {
        // Skip node if it is the dummy node created in BVH (nodeIndex == 1), 
        // or the level of this node is not the required one.
        if (nodeIndex == 1 || levels[nodeIndex] != static_cast<uint32_t>(level)) {
            // This is the dummy node created in BVH constructor.
            continue;
        }

        drawAABB(m_nodes[nodeIndex].aabb, DrawMode::Wireframe, color, transparency);
    }
}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    const glm::vec3 aabbColor { 0.5f, 0.5f, 1.0f};
    const float aabbTransparency = 0.8f;

    // Each primitive in the leaf will have a slightly different color.
    // This makes it easier to see the rectangles themselves.
    glm::vec3 triangleColor { 0.5f, 1.0f, 0.4f };
    const glm::vec3 triangleColorIncrement { 0.0f, -0.2f, 0.2f };

    int currentLeafIndex = 0;
    for (int i = 0; i < m_nodes.size(); i++) {
        // Skip if it is not a leaf, or it is a dummy node (added in BVH constructor).
        if (i == 1 || !m_nodes[i].isLeaf()) {
            continue;
        }

        // m_nodes[i] is a leaf.
        currentLeafIndex++;

        // Skip if it is not the required leaf.
        if (leafIndex != currentLeafIndex) {
            continue;
        }

        // This node (which is a leaf) should be drawn.
        // Firstly, draw AABB of this leaf.
        drawAABB(m_nodes[i].aabb, DrawMode::Wireframe, aabbColor, aabbTransparency);

        // Draw all primitives contained in this leaf.
        for (uint32_t primitiveId = m_nodes[i].primitiveOffset();
             primitiveId < m_nodes[i].primitiveOffset() + m_nodes[i].primitiveCount();
             primitiveId++) {
            Primitive primitive = m_primitives[primitiveId];
            glColor3f(triangleColor[0], triangleColor[1], triangleColor[2]);
            drawTriangle(primitive.v0, primitive.v1, primitive.v2);
            triangleColor += triangleColorIncrement;
        }

        break; // No need to continue the loop.
    }
}