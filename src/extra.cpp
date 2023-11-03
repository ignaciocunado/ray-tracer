#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "texture.h"
#include "draw.h"
#include <framework/trackball.h>

// Helper method for renderImageWithDepthOfField(...).
// Calculates the position of the point which is on the focal plane, and which is hit by cameraRay.
glm::vec3 getPointOfFocus(const Trackball& camera, const Ray& cameraRay, const float focalDistance)
{
    const float cosAngle = glm::dot(cameraRay.direction, camera.forward()); // cosine of the angle between these two vectors

    // The following line is get due to: cosAngle = focalDistance / curDistance
    const float curDistance = focalDistance / cosAngle; // distance for the cameraRay to reach focal plane

    return cameraRay.origin + curDistance * cameraRay.direction;
}

// Helper method for renderImageWithDepthOfField(...).
// Calculates rays for specified pixel (x, y) and renders.
void renderImagePixelWithDepthOfField(RenderState& state, const Trackball& camera, Screen& screen, const glm::ivec2& pixel)
{
    float focalDistance = state.features.extra.depthOfFieldDistance;
    float squareLength = state.features.extra.depthOfFieldSquareLength;
    uint32_t numOfSamples = state.features.extra.numDepthOfFieldSamples;

    // Generate rays from camera for this pixel (x, y). (Similarly to renderImage(...))
    std::vector<Ray> generatedRays = generatePixelRays(state, camera, pixel, screen.resolution());

    std::vector<Ray> finalRays;
    finalRays.reserve(numOfSamples * generatedRays.size());

    // For each of the generated camera ray new rays will be created for the Depth of field effect.
    for (const Ray& cameraRay : generatedRays) {
        const glm::vec3 pointOfFocus = getPointOfFocus(camera, cameraRay, focalDistance);

        // Generate many new rays.
        for (uint32_t i = 0; i < numOfSamples; i++) {
            // Generate new ray with origin slightly moved, but the direction still
            // directed towards point of focus.

            const glm::vec2 randomOffset = (state.sampler.next_2d() - 0.5f) * squareLength;
            // randomOffset is added to the camera plane which is orthogonal to camera.forward().
            // The orthogonal basis (two orthogonal unit vectors) can be easily retrieved by taking camera.up() and camera.left().
            const glm::vec3 newOrigin = cameraRay.origin + randomOffset[0] * camera.up() + randomOffset[1] * camera.left();
            const glm::vec3 newDirection = glm::normalize(pointOfFocus - newOrigin);

            finalRays.push_back(Ray { .origin = newOrigin, .direction = newDirection, .t = cameraRay.t });
        }
    }

    auto L = renderRays(state, finalRays);
    screen.setPixel(pixel.x, pixel.y, L);
}

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

#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
    // Loop through each pixel.
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            RenderState state = {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
            };

            renderImagePixelWithDepthOfField(state, camera, screen, { x, y });
        }
    }
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
}

// Helper function for factorial
long long factorial(int n)
{
    long long f = 1;

    for (int i = 1; i <= n; i++) {
        f *= i;
    }

    return f;
}

void computeGaussianFilter(std::vector<std::vector<float>>& filter, int k)
{
    // Calculate the horizontal values
    for (int y = 0; y < k; y++) {
        long double total = 0;

        for (int x = 0; x < k; x++) {
            float val = float(factorial(k) / (factorial(x) * factorial(k - x)));
            total += val;
            filter[x][y] = val;
        }

        for (int x = 0; x < k; x++) {
            filter[x][y] /= total;
        }
    }

    // Calculate the vertical values
    for (int x = 0; x < k; x++) {
        long long total = 0;

        for (int y = 0; y < k; y++) {
            int val = int(factorial(k) / (factorial(y) * factorial(k - y)));
            total += val;
            filter[x][y] = val;
        }

        for (int y = 0; y < k; y++) {
            filter[x][y] /= total;
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

    glm::ivec2 resolution = image.resolution();
    int width = resolution.x;
    int height = resolution.y;
    int k = int(features.extra.bloomFilterSize);
    std::vector<std::vector<float>> filter(k, std::vector<float>(k));
    std::vector<glm::vec3> originalPixels = image.pixels();
    ;

    // Compute the gaussian filter
    computeGaussianFilter(filter, k);

    // Apply the k x k filter to the image, leave the border pixels alone
    for (int x = k - 2; x < width - k + 2; x++) {
        for (int y = k - 2; y < height - k + 2; y++) {
            glm::vec3 color = glm::vec3(0.f);

            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    color += filter[i][j] * originalPixels[image.indexAt(x + i, y + j)];
                }
            }

            image.setPixel(x, y, color);
        }
    }
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
    glm::vec3 glossyColor = glm::vec3 { 0.f };
    // Generate an initial specular ray, and base secondary glossies on this ray
    glm::vec3 intersection = ray.origin + ray.direction * ray.t;
    Ray initialSpecularRay = generateReflectionRay(ray, hitInfo);

    // Construct the orthonormal basis for the specular ray
    glm::vec3 arbitraryVector = glm::vec3(1, 0, 0);
    if (glm::length(arbitraryVector - initialSpecularRay.direction) < 0.01f) {
        arbitraryVector = glm::vec3(0, 1, 0);
    }

    glm::vec3 u = glm::normalize(glm::cross(arbitraryVector, initialSpecularRay.direction));
    glm::vec3 v = glm::normalize(glm::cross(initialSpecularRay.direction, u));

    int numSamples = (int)state.features.extra.numGlossySamples;
    float diskRadius = state.features.extra.glossyExponent * hitInfo.material.shininess / 64;

    for (int i = 0; i < numSamples; i++) {
        // Generate a sampled point on the disk
        glm::vec2 samples = state.sampler.next_2d();
        float sampledDiskRadius = diskRadius * samples.x;
        float theta = glm::two_pi<float>() * samples.y;

        // Get the sampled coordinates
        glm::vec3 sampledU = u * sampledDiskRadius * cos(theta);
        glm::vec3 sampledV = v * sampledDiskRadius * sin(theta);

        glm::vec3 sampledDirection = glm::normalize(initialSpecularRay.direction + sampledU + sampledV);

        // Generate a ray and render it
        Ray glossyRay = Ray { .origin = intersection + sampledDirection * 1.0e-5f, .direction = sampledDirection, .t = std::numeric_limits<float>::max() };
        drawRay(glossyRay, glm::vec3 { 0, 1, 1 });
        glossyColor += renderRay(state, glossyRay, rayDepth + 1);
    }

    // Add the color to the hitColor
    hitColor += (glossyColor / float(numSamples)) * hitInfo.material.ks;
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
    if (!state.features.extra.enableEnvironmentMap) {
        return glm::vec3(0.0f); // Black color.
    }

    // Index of the face (image) of the cube (of environment map textures). Not yet known.
    // Indices meaning (order) is listed in scene.h file.
    int chosenImageId = 0;
    // Texture coordinates. Not yet known.
    float u = 0.0f;
    float v = 0.0f;

    // Calculate which side of the cube the ray is facing.
    const float x = ray.direction.x;
    const float y = ray.direction.y;
    const float z = ray.direction.z;
    const float absX = abs(x);
    const float absY = abs(y);
    const float absZ = abs(z);

    // The side of the cube can be figured out by checking which coordinate has biggest value, and checking its sign.
    if (absX >= absY && absX >= absZ) {
        // X is biggest, so either positive or negative X.
        if (x > 0.0f) {
            // Facing positive x direction.
            chosenImageId = 0;
            u = z;
            v = y;
        } else {
            // Facing negative x direction.
            chosenImageId = 1;
            u = -z;
            v = y;
        }
    } else if (absY >= absZ) {
        // Y is biggest, so either positive or negative Y.
        if (y > 0.0f) {
            // Facing positive y direction.
            chosenImageId = 2;
            u = -x;
            v = -z;
        } else {
            // Facing negative y direction.
            chosenImageId = 3;
            u = -x;
            v = z;
        }
    } else {
        // Z is biggest, so either positive or negative Z.
        if (z > 0.0f) {
            // Facing positive z direction.
            chosenImageId = 4;
            u = -x;
            v = y;
        } else {
            // Facing negative z direction.
            chosenImageId = 5;
            u = x;
            v = y;
        }
    }

    // Check if texture is present.
    if (!state.scene.environmentMapTextures[chosenImageId]) {
        return glm::vec3(0.0f); // Black color.
    }

    // Normalize (u, v) coordinates to be between [0; 1]x[0; 1].
    const float biggestAbsCoord = std::max(absX, std::max(absY, absZ));
    u = (u / biggestAbsCoord + 1.0f) * 0.5f;
    v = (v / biggestAbsCoord + 1.0f) * 0.5f;

    const Image& chosenImage = *(state.scene.environmentMapTextures[chosenImageId]);
    const glm::vec2 texCoords(u, v);

    if (state.features.enableBilinearTextureFiltering) {
        return sampleTextureBilinear(chosenImage, texCoords);
    } else {
        return sampleTextureNearest(chosenImage, texCoords);
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