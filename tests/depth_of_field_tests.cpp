// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "shading.h"
#include "extra.h"
#include <limits>
#include <framework/trackball.h>
#include <framework/window.h>

#include "intersect.h"

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()

namespace dof_tests {
void checkVec3(const glm::vec3& actual, const glm::vec3& expected, const float eps = 1e-9f)
{
    for (int coord = 0; coord < 3; coord++) {
        REQUIRE_THAT(actual[coord], Catch::Matchers::WithinAbs(expected[coord], eps));
    }
}
}

TEST_CASE("Calculating Focus points") {
    const glm::ivec2 windowSize { 800, 800 };
    Window window { "Window", windowSize, OpenGLVersion::GL2, false };
    Screen screen { windowSize, false };
    Trackball camera { &window, glm::radians(50.0f), 3.0f};
    
    // Use GeoGebra 3D calculator to plot the values.
    // https://www.geogebra.org/m/d7yjyetq (do NOT share this link to anyone!)
    const glm::vec3 lookAt(1, 2, 5);
    const glm::vec3 rotationsInRad(1.0f, 1.5f, 2.0f);
    const float distanceFromLookAt = 3.0f;
    camera.setCamera(lookAt, rotationsInRad, distanceFromLookAt);
    const float focalLength = 4.0f;

    const float epsBig = 1e-6f; // quite a big epsilon due to inaccuracies of converting between degrees and radians (inaccuracy of PI constant)

    {
        // Check that camera is set up as expected.
        dof_tests::checkVec3(camera.position(), glm::vec3(-0.622596637f, -0.520720849f, 4.885341580f), epsBig);
        dof_tests::checkVec3(camera.forward(), glm::vec3(0.540865545f, 0.840240283f, 0.038219473f), epsBig);
        dof_tests::checkVec3(camera.left(), glm::vec3(-0.029437062f, 0.064321155f, -0.99749486), epsBig);
        dof_tests::checkVec3(camera.up(), glm::vec3(-0.840593790f, 0.538385601f, 0.059523302f), epsBig);
    }

    SECTION("Focal point of the center")
    {
        const Ray ray = { .origin = camera.position(), .direction = camera.forward() };
        const glm::vec3 result = getPointOfFocus(camera, ray, focalLength);
        dof_tests::checkVec3(result, glm::vec3(1.540865545f, 2.840240283f, 5.038219473f), epsBig);
        
    }

    SECTION("Some focal point #1 - geogebra B")
    {
        const Ray ray = { .origin = camera.position(), .direction = glm::vec3(0.707073733f, 0.675492054f, -0.209182263f) };
        const glm::vec3 result = getPointOfFocus(camera, ray, focalLength);
        dof_tests::checkVec3(result, glm::vec3(2.379799376f, 2.347572165f, 3.997106076f), epsBig);
    }

    SECTION("Some focal point #2 - geogebra C")
    {
        const Ray ray = { .origin = camera.position(), .direction = glm::vec3(0.007410816f, 0.919139747f, -0.393861909f) };
        const glm::vec3 result = getPointOfFocus(camera, ray, focalLength);
        dof_tests::checkVec3(result, glm::vec3(-0.583656555f, 4.308892342f, 2.815797010f), epsBig);
    }

    SECTION("Some focal point #3 - geogebra D")
    {
        const Ray ray = { .origin = camera.position(), .direction = glm::vec3(0.302661460f, 0.791110773f, 0.531544714f) };
        const glm::vec3 result = getPointOfFocus(camera, ray, focalLength);
        dof_tests::checkVec3(result, glm::vec3(0.803810963f, 3.207690485f, 7.390448834f), epsBig);
    }
}
