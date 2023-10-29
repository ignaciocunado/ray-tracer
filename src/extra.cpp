#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

// Helper method for renderImageWithDepthOfField(...).
// Calculates the position of the point which is on the focal plane, and which is hit by cameraRay.
glm::vec3 getPointOfFocus(const Trackball& camera, const Ray& cameraRay, const float focalDistance) {
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

            renderImagePixelWithDepthOfField(state, camera, screen, {x, y});
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