// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "light.h"
#include <limits>

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()

// In this file you can add your own unit tests using the Catch2 library.
// You can find the documentation of Catch2 at the following link:
// https://github.com/catchorg/Catch2/blob/devel/docs/assertions.md
//
// These tests are only to help you verify that your code is correct.
// You don't have to hand them in; we will not consider them when grading.
//

TEST_CASE("Sample Segment")
{
    float epsilon = 0.0001f;

    struct TestCase {
        float sample;
        SegmentLight light;
        glm::vec3 expectedPosition;
        glm::vec3 expectedColor;
    };

    TestCase testCases[] = {
        {0.0f, {glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0.2f, 0.2f, 0.2f), glm::vec3(0.8f, 0.8f, 0.8f)},
            glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.2f, 0.2f, 0.2f)},
        {0.3f, {glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0.2f, 0.2f, 0.2f), glm::vec3(0.8f, 0.8f, 0.8f)},
            glm::vec3(0.3f, 0.3f, 0.3f), glm::vec3(0.38f, 0.38f, 0.38f)},
        {0.5f, {glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0.2f, 0.2f, 0.2f), glm::vec3(0.8f, 0.8f, 0.8f)},
            glm::vec3(0.5f, 0.5f, 0.5f), glm::vec3(0.5f, 0.5f, 0.5f)},
        {1.0f, {glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0.2f, 0.2f, 0.2f), glm::vec3(0.8f, 0.8f, 0.8f)},
            glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0.8f, 0.8f, 0.8f)},
    };

    for (const auto& testCase : testCases) {
        glm::vec3 position, color;
        sampleSegmentLight(testCase.sample, testCase.light, position, color);

        // Check if the actual results match the expected results
        REQUIRE(glm::all(glm::epsilonEqual(position, testCase.expectedPosition, epsilon)));
        REQUIRE(glm::all(glm::epsilonEqual(color, testCase.expectedColor, epsilon)));
    }
}

TEST_CASE("Barycentric Coordinates Calculation") {
    // Define the edges e1 and e2
    glm::vec3 e1(1.0f, 0.0f, 0.0f);
    glm::vec3 e2(0.0f, 1.0f, 0.0f);
    
    SECTION("Test at sample (0.5, 0.5)") {
        glm::vec2 sample(0.5f, 0.5f);
        glm::vec4 barycentric = calculateBarycentricCoordinates(e1, e2, sample);
        REQUIRE(barycentric.x == Catch::Approx(0.25f));
        REQUIRE(barycentric.y == Catch::Approx(0.25f));
        REQUIRE(barycentric.z == Catch::Approx(0.25f));
        REQUIRE(barycentric.w == Catch::Approx(0.25f));
    }

    SECTION("Test at sample (0.0, 0.0)") {
        glm::vec2 sample(0.0f, 0.0f);
        glm::vec4 barycentric = calculateBarycentricCoordinates(e1, e2, sample);
        REQUIRE(barycentric.x == Catch::Approx(0.0f));
        REQUIRE(barycentric.y == Catch::Approx(0.0f));
        REQUIRE(barycentric.z == Catch::Approx(0.0f));
        REQUIRE(barycentric.w == Catch::Approx(1.0f));
    }

    SECTION("Test at sample (1.0, 0.0)") {
        glm::vec2 sample(1.0f, 0.0f);
        glm::vec4 barycentric = calculateBarycentricCoordinates(e1, e2, sample);
        REQUIRE(barycentric.x == Catch::Approx(0.0f));
        REQUIRE(barycentric.y == Catch::Approx(0.0f));
        REQUIRE(barycentric.z == Catch::Approx(1.0f));
        REQUIRE(barycentric.w == Catch::Approx(0.0f));
    }

    SECTION("Test at sample (0.0, 1.0)") {
        glm::vec2 sample(0.0f, 1.0f);
        glm::vec4 barycentric = calculateBarycentricCoordinates(e1, e2, sample);
        REQUIRE(barycentric.x == Catch::Approx(0.0f));
        REQUIRE(barycentric.y == Catch::Approx(1.0f));
        REQUIRE(barycentric.z == Catch::Approx(0.0f));
        REQUIRE(barycentric.w == Catch::Approx(0.0f));
    }

    SECTION("Test at sample (1.0, 1.0)") {
        glm::vec2 sample(1.0f, 1.0f);
        glm::vec4 barycentric = calculateBarycentricCoordinates(e1, e2, sample);
        REQUIRE(barycentric.x == Catch::Approx(1.0f));
        REQUIRE(barycentric.y == Catch::Approx(0.0f));
        REQUIRE(barycentric.z == Catch::Approx(0.0f));
        REQUIRE(barycentric.w == Catch::Approx(0.0f));
    }
}

TEST_CASE("sampleParallelogramLight function") {
    // Define a ParallelogramLight instance for testing
    ParallelogramLight light;
    light.v0 = glm::vec3(0.0f, 0.0f, 0.0f);
    light.edge01 = glm::vec3(1.0f, 0.0f, 0.0f);
    light.edge02 = glm::vec3(0.0f, 1.0f, 0.0f);
    light.color0 = glm::vec3(1.0f, 0.0f, 0.0f);
    light.color1 = glm::vec3(0.0f, 1.0f, 0.0f);
    light.color2 = glm::vec3(0.0f, 0.0f, 1.0f);
    light.color3 = glm::vec3(1.0f, 1.0f, 1.0f);

    SECTION("Test at sample (0.0, 0.0)") {
        glm::vec2 sample(0.0f, 0.0f);
        glm::vec3 position, color;
        sampleParallelogramLight(sample, light, position, color);

        REQUIRE(position == glm::vec3(0.0f, 0.0f, 0.0f));
        REQUIRE(color == glm::vec3(1.0f, 0.0f, 0.0f));
    }

    SECTION("Test at sample (0.3, 0.7)") {
        glm::vec2 sample(0.3f, 0.7f);
        glm::vec3 position, color;
        sampleParallelogramLight(sample, light, position, color);

        glm::vec3 expectedPosition(0.3f, 0.7f, 0.0f);
        glm::vec3 expectedColor = glm::mix(glm::mix(light.color0, light.color1, 0.3f), glm::mix(light.color2, light.color3, 0.3f), 0.7f);

        REQUIRE(position == expectedPosition);
        REQUIRE(color == expectedColor);
    }

    SECTION("Test at sample (0.5, 0.5)") {
        glm::vec2 sample(0.5f, 0.5f);
        glm::vec3 position, color;
        sampleParallelogramLight(sample, light, position, color);

        REQUIRE(position == glm::vec3(0.5f, 0.5f, 0.0f));
        REQUIRE(color == glm::vec3(0.5f, 0.5f, 0.5f));
    }

    SECTION("Test at sample (1.0, 1.0)") {
        glm::vec2 sample(1.0f, 1.0f);
        glm::vec3 position, color;
        sampleParallelogramLight(sample, light, position, color);

        REQUIRE(position == glm::vec3(1.0f, 1.0f, 0.0f));
        REQUIRE(color == glm::vec3(1.0f, 1.0f, 1.0f));
    }
}