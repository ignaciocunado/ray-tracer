// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "shading.h"
#include "interpolate.h"
#include <limits>

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()

void vectorCheck(const glm::vec3 result, const glm::vec3 expected, const float eps = 1e-6)
{
    CHECK_THAT(result[0], Catch::Matchers::WithinAbs(expected[0], eps));
    CHECK_THAT(result[1], Catch::Matchers::WithinAbs(expected[1], eps));
    CHECK_THAT(result[2], Catch::Matchers::WithinAbs(expected[2], eps));
}

TEST_CASE("barycentric") {
    glm::vec3 v0 { 3, 4, 5 };
    glm::vec3 v1 { 0, 5, 3 };
    glm::vec3 v2 { 3, 0, -2 };
    glm::vec3 p { 2, 3, 2 };
    glm::vec3 expected { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f };
    vectorCheck(computeBarycentricCoord(v0, v1, v2, p), expected);
}

TEST_CASE("barycentricEdge")
{
    glm::vec3 v0 { 0, 0, 0 };
    glm::vec3 v1 { 1, 0, 0 };
    glm::vec3 v2 { 0.4, 1, 0 };
    glm::vec3 p { 0.6, 0, 0 };
    glm::vec3 expected { 0.4f, 0.6f, 0.0f };
    vectorCheck(computeBarycentricCoord(v0, v1, v2, p), expected);
}

TEST_CASE("barycentricVertex")
{
    glm::vec3 v0 { 0, 0, 0 };
    glm::vec3 v1 { 1, 0, 0 };
    glm::vec3 v2 { 0.4, 1, 0 };
    glm::vec3 p { 0.4, 1, 0 };
    glm::vec3 expected { 0.0f, 0.0f, 1.0f };
    vectorCheck(computeBarycentricCoord(v0, v1, v2, p), expected);
}

TEST_CASE("interpolateNormals")
{
    glm::vec3 n0 { 1,0,0 };
    glm::vec3 n1 { 0, 1, 0 };
    glm::vec3 n2 { 0, 0, 1 };
    glm::vec3 bar { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f };
    glm::vec3 expected { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f };
    vectorCheck(interpolateNormal(n0, n1, n2, bar), expected);
}
