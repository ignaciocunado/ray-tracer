#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded linear gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            { 0.1f, glm::vec3(215.f / 256.f, 210.f / 256.f, 203.f / 256.f) },
            { 0.22f, glm::vec3(250.f / 256.f, 250.f / 256.f, 240.f / 256.f) },
            { 0.5f, glm::vec3(145.f / 256.f, 170.f / 256.f, 175.f / 256.f) },
            { 0.78f, glm::vec3(255.f / 256.f, 250.f / 256.f, 205.f / 256.f) },
            { 0.9f, glm::vec3(170.f / 256.f, 170.f / 256.f, 170.f / 256.f) },
        }
    };

    if (state.features.enableShading) {
        switch (state.features.shadingModel) {
            case ShadingModel::Lambertian:
                return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::Phong:
                return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::BlinnPhong:
                return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo);
            case ShadingModel::LinearGradient:
                return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
        };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a Lambertian diffuse shading, returning the reflected light towards the target.
glm::vec3 computeLambertianModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Compute the dot product between the normal and the light direction
    const float dot = glm::dot(hitInfo.normal, lightDirection);

    // Check if the diffuse component should be computed
    if (dot < 0) {
        return glm::vec3(0.f);
    }

    const glm::vec3 materialKd = sampleMaterialKd(state, hitInfo);

    // Compute the diffuse component
    return lightColor * materialKd * dot;
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Compute the diffuse component
    const glm::vec3 diffuse = computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);

    // Compute the reflected ray direction and the dot product between the reflected ray direction and the camera direction
    const glm::vec3 reflected = glm::reflect(lightDirection, hitInfo.normal);
    const float dot = glm::dot(reflected, cameraDirection);

    // Check if the specular component should be computed
    if (dot < 0) {
        return diffuse;
    }

    // Compute the specular component
    const glm::vec3 materialKs = hitInfo.material.ks;
    const float shininess = hitInfo.material.shininess;
    const glm::vec3 specular = lightColor * materialKs * glm::pow(dot, shininess);

    // Return the sum of the diffuse and specular components (ambient component is ignored)
    return diffuse + specular;
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the Blinn-Phong Model returning the reflected light towards the target.
// Note: materials do not have an ambient component, so you can ignore this.
// Note: use `sampleMaterialKd` instead of material.kd to automatically forward to texture
//       sampling if a material texture is available!
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - return;          the result of shading along the cameraDirection vector
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBlinnPhongModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Compute the diffuse component
    const glm::vec3 diffuse = computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo);

    // Check if the specular component should be computed
    if (glm::dot(hitInfo.normal, lightDirection) < 0) {
        return diffuse;
    }

    // Compute the specular component
    const glm::vec3 H = glm::normalize(lightDirection + cameraDirection);
    const float dot = glm::dot(hitInfo.normal, H);
    const float shininess = hitInfo.material.shininess;
    const glm::vec3 materialKs = hitInfo.material.ks;
    const glm::vec3 specular = lightColor * materialKs * glm::pow(dot, shininess);

    return diffuse + specular;
}

// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    // Check if the gradient is empty
    if (components.empty()) {
        return glm::vec3(0.f);
    }

    // Handle special cases when ti is outside the gradient's range
    if (ti <= components.front().t) {
        return components.front().color;
    } else if (ti >= components.back().t) {
        return components.back().color;
    }

    // Find the two components between which ti lies
    auto it = std::lower_bound(components.begin(), components.end(), ti,[](const Component& component, float value) {
        return component.t <= value;
    });

    const Component& component1 = *(it - 1);
    const Component& component2 = *it;

    // Compute the interpolation factor
    const float t = (ti - component1.t) / (component2.t - component1.t);

    // Interpolate between the two components
    return glm::mix(component1.color, component2.color, t);
}

// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    // Compute the cosine of theta and sample the gradient
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    glm::vec3 materialKd = gradient.sample(cos_theta);

    // Calculate the diffuse component
    return lightColor * materialKd;
}