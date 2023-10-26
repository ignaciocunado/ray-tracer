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
    /*static glm::mat4 spliceMat(float t)
    {
        if (t >= 1) {
            t = 0;
        }
        glm::vec3 p0 { 0, 0, 0 };
        glm::vec3 p1 { 0, 5, 0 };
        glm::vec3 p2 { 5, 5, 0 };
        glm::vec3 p3 { 5, 0, 0 };
        glm::vec3 posBezier = ((1.0f - t) * (1.0f - t) * (1.0f - t) * p0) + (3 * (1.0f - t) * (1.0f - t) * t * p1) + (3 * (1.0f - t) * t * t * p2) * (t * t * t * p3);
        return glm::translate(glm::identity<glm::mat4>(), posBezier);
    }*/

    // Either directly render the image, or pass through to extra.h methods
    if (features.extra.enableDepthOfField) {
        renderImageWithDepthOfField(scene, bvh, features, camera, screen);
    } else if (features.extra.enableMotionBlur) {
        renderImageWithMotionBlur(scene, bvh, features, camera, screen);
    } else {
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
                auto rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
                auto L = renderRays(state, rays);
                screen.setPixel(x, y, L);
            }
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