#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

    int samples = features.extra.motionBlurSamples;

    #ifdef NDEBUG // Enable multi threading in Release mode
    #pragma omp parallel for schedule(guided)
    #endif
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
            // Note; we seed the sampler for consistenct behavior across frames
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };
            
            glm::vec3 L;
            for (int i = 0; i < samples; i++) {
                float time = state.sampler.next_1d();

                std::vector<Sphere> newSpheres;
                for (int j = 0; j < scene.spheres.size(); j++) {
                    Sphere change = scene.spheres[j];
                    glm::vec3 center = change.center;
                    glm::mat4 transform = spliceMat(time, center);
                    glm::vec4 newCenter = transform * glm::vec4 { center[0], center[1], center[2], 1 };
                    glm::vec3 newCenter3 = { newCenter[0] / newCenter[3], newCenter[1] / newCenter[3], newCenter[2] / newCenter[3] };
                    Sphere newSphere = {
                        .center = newCenter3,
                        .radius = change.radius,
                        .material = change.material
                    };
                    newSpheres.push_back(newSphere);
                }

                std::vector<Mesh> newMeshes;
                for (int j = 0; j < scene.meshes.size(); j++) {
                    Mesh change = scene.meshes[j];
                    std::vector<Vertex> vertices = change.vertices;
                    std::vector<Vertex> newVertices;
                    for (int k = 0; k < vertices.size(); k++) {
                        Vertex v = vertices[k];
                        glm::vec3 pos = v.position;
                        glm::mat4 transform = spliceMat(time, pos);
                        glm::vec4 newPos = transform * glm::vec4 { pos[0], pos[1], pos[2], 1 };
                        glm::vec3 newPos3 = { newPos[0] / newPos[3], newPos[1] / newPos[3], newPos[2] / newPos[3] };
                        Vertex newVertex = {
                            .position = newPos3,
                            .normal = v.normal,
                            .texCoord = v.texCoord
                        };
                        newVertices.push_back(newVertex);
                    }
                    Mesh newMesh = {
                        .vertices = newVertices,
                        .triangles = change.triangles,
                        .material = change.material
                    };
                    newMeshes.push_back(newMesh);
                }

                Scene newScene
                {
                    .type = scene.type,
                    .meshes = newMeshes,
                    .spheres = newSpheres,
                    .lights = scene.lights

                };
                BVH bvh = BVH(newScene, features);
                RenderState newState = {
                    .scene = newScene,
                    .features = state.features,
                    .bvh = bvh,
                    .sampler = state.sampler
                };

                auto rays = generatePixelRays(newState, camera, { x, y }, screen.resolution());
                L += renderRays(newState, rays);
            }
            L = L / (float) samples;
            screen.setPixel(x, y, L);
        }
    }
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    return 0; // This is clearly not the solution
}

// This method calculates the position of a point along a 4th degree bezier curve at time T
glm::mat4 spliceMat(float t, glm::vec3 currentCenter)
{
    glm::vec3 p0 = glm::vec3(0, 0, 0) + currentCenter;
    glm::vec3 p1 = glm::vec3(0, 5, 2) + currentCenter;
    glm::vec3 p2 = glm::vec3(5, 5, -2) + currentCenter;
    glm::vec3 p3 = glm::vec3(5, 0, 0) + currentCenter;
    glm::vec3 p4 = glm::vec3(2.5, 10, 4) + currentCenter;

    float oneMinusT = 1.0f - t;
    float oneMinusTSquared = oneMinusT * oneMinusT;
    float tSquared = t * t;
    float tCubed = tSquared * t;

    glm::vec3 posBezier = (oneMinusTSquared * oneMinusT * oneMinusT * p0) + (4.0f * oneMinusTSquared * oneMinusT * t * p1) + (6.0f * oneMinusTSquared * tSquared * p2) + (4.0f * oneMinusT * tCubed * p3) + (tSquared * tCubed * p4);

    return glm::translate(glm::mat4(1.0f), posBezier);
}