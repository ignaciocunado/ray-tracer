#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
#include "recursive.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    glm::vec3 p0 = light.endpoint0;
    glm::vec3 p1 = light.endpoint1;

    glm::vec3 c0 = light.color0;
    glm::vec3 c1 = light.color1;

    position = glm::mix(p0, p1, sample);
    color = glm::mix(c0, c1, sample);
}

glm::vec4 calculateBarycentricCoordinates(const glm::vec3& e1, const glm::vec3& e2, const glm::vec2& sample) {
    float x = sample.x, y = sample.y;

    float area0 = glm::length(glm::cross(e1 * x, e2 * y));
    float area1 = glm::length(glm::cross(e1 * (1 - x), e2 * y));
    float area2 = glm::length(glm::cross(e1 * x, e2 * (1 - y)));
    float area3 = glm::length(glm::cross(e1 * (1 - x), e2 * (1 - y)));

    float totalArea = area0 + area1 + area2 + area3;

    return glm::vec4(area0 / totalArea, area1 / totalArea, area2 / totalArea, area3 / totalArea);
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    glm::vec3 p0 = light.v0;
    glm::vec3 p1 = light.v0 + light.edge01;
    glm::vec3 p2 = light.v0 + light.edge02;
    position = p0 + sample.x * light.edge01 + sample.y * light.edge02;

    glm::vec3 c0 = light.color0;
    glm::vec3 c1 = light.color1;
    glm::vec3 c2 = light.color2;
    glm::vec3 c3 = light.color3;

    glm::vec4 barycentricCoordinates = calculateBarycentricCoordinates(light.edge01, light.edge02, sample);
    color = barycentricCoordinates.w * c0 + barycentricCoordinates.z * c1 + barycentricCoordinates.y * c2 + barycentricCoordinates.x * c3;
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Create a ray from the light position to the hit point
        glm::vec3 hitPoint = ray.origin + ray.t * ray.direction;
        glm::vec3 shadowRayDirection = glm::normalize(hitPoint - lightPosition);
        float shadowRayLength = glm::length(hitPoint - lightPosition);
        Ray shadowRay = {lightPosition, shadowRayDirection, shadowRayLength};
        HitInfo shadowHitInfo;

        // Use the BVH intersect function to check for occlusions
        if (state.bvh.intersect(state, shadowRay, shadowHitInfo)) {
            // Check if the shadow ray hit an object before the light source
            return abs(shadowRay.t - shadowRayLength) <= 0.0001f;
        }

        // The shadow ray didn't hit any objects, so the light is visible
        return true;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 hitPoint = ray.origin + ray.t * ray.direction;
    glm::vec3 shadowRayDirection = glm::normalize(hitPoint - lightPosition);
    float shadowRayLength = glm::length(hitPoint - lightPosition);
    Ray shadowRay = {lightPosition, shadowRayDirection, shadowRayLength};
    HitInfo shadowHitInfo;
    glm::vec3 finalColor {0.0f};
    float transparency = 1.0f;
    Material currentMaterial = hitInfo.material;

    // Check if the shadow ray hit an object before the hitPoint
    while (state.bvh.intersect(state, shadowRay, shadowHitInfo) && !glm::all(glm::epsilonEqual(shadowRay.origin + shadowRay.t * shadowRay.direction, hitPoint, 0.0001f)) && shadowHitInfo.material.transparency != 1.0f) {
        // Ray hit a transparent object, generate a pass-through ray and save the current material
        currentMaterial = shadowHitInfo.material;
        transparency *= currentMaterial.transparency;
        shadowRay = generatePassthroughRay(shadowRay, shadowHitInfo);
    }

    // Check if the ray reached the object
    if (glm::all(glm::epsilonEqual(shadowRay.origin + shadowRay.t * shadowRay.direction, hitPoint, 0.0001f))) {
        finalColor = lightColor * currentMaterial.kd * transparency;
    }

    return finalColor;
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 visibleColor = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);

    // If the light is not visible, return black
    if (visibleColor == glm::vec3(0)) {
        return glm::vec3(0);
    }

    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;
    HitInfo lightHitInfo = hitInfo;

    if (state.features.enableTransparency && hitInfo.material.transparency < 1.0f && glm::dot(hitInfo.normal, l) < 0.0f) {
        lightHitInfo.normal = -hitInfo.normal;
    }

    // Compute the shading for the point light using the visible color
    return computeShading(state, v, l, visibleColor, lightHitInfo);
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    glm::vec3 Lo { 0.0f };

    for (uint32_t i = 0; i < numSamples; i++) {
        // - sample the segment light
        glm::vec3 position, color;
        sampleSegmentLight(state.sampler.next_1d(), light, position, color);

        Lo += computeContributionPointLight(state, PointLight {position, color}, ray, hitInfo);
    }

    return Lo * (1.0f / (float)numSamples);
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    glm::vec3 Lo { 0.0f };

    for (uint32_t i = 0; i < numSamples; i++) {
        // - sample the parallelogram light
        glm::vec3 position, color;
        sampleParallelogramLight(state.sampler.next_2d(), light, position, color);

        Lo += computeContributionPointLight(state, PointLight {position, color}, ray, hitInfo);
    }

    return Lo * (1.0f / (float)numSamples);
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }

    return glm::clamp(Lo, 0.0f, 1.0f);
}