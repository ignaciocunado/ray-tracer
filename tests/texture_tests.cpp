// Put your includes here
#include "texture.h"
#include <limits>

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()

void vec3Check(const glm::vec3 result, const glm::vec3 expected, const float eps = 1e-6)
{
    CHECK_THAT(result[0], Catch::Matchers::WithinAbs(expected[0], eps));
    CHECK_THAT(result[1], Catch::Matchers::WithinAbs(expected[1], eps));
    CHECK_THAT(result[2], Catch::Matchers::WithinAbs(expected[2], eps));
}

TEST_CASE("nearest1") {
    Image image = Image("../../../../tests/image.png");
    image.height = 3;
    image.width = 3;
    image.pixels = std::vector<glm::vec3> {
        glm::vec3 { 1 }, glm::vec3 { 2 }, glm::vec3 { 3 }, glm::vec3 { 4 }, glm::vec3 { 5 }, glm::vec3 { 6 }, glm::vec3 { 7 }, glm::vec3 { 8 }, glm::vec3 { 9 }
    };
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0 }), glm::vec3 { 7 });
}

TEST_CASE("nearest2")
{
    Image image = Image("../../../../tests/image.png");
    image.height = 3;
    image.width = 3;
    image.pixels = std::vector<glm::vec3> {
        glm::vec3 { 1 }, glm::vec3 { 2 }, glm::vec3 { 3 }, glm::vec3 { 4 }, glm::vec3 { 5 }, glm::vec3 { 6 }, glm::vec3 { 7 }, glm::vec3 { 8 }, glm::vec3 { 9 }
    };
    vec3Check(sampleTextureNearest(image, glm::vec2 { 1 ,0 }), glm::vec3 { 9 });
}

TEST_CASE("nearest3")
{
    Image image = Image("../../../../tests/image.png");
    image.height = 3;
    image.width = 3;
    image.pixels = std::vector<glm::vec3> {
        glm::vec3 { 1 }, glm::vec3 { 2 }, glm::vec3 { 3 }, glm::vec3 { 4 }, glm::vec3 { 5 }, glm::vec3 { 6 }, glm::vec3 { 7 }, glm::vec3 { 8 }, glm::vec3 { 9 }
    };
    vec3Check(sampleTextureNearest(image, glm::vec2 { 1, 1 }), glm::vec3 { 3 });
}

TEST_CASE("nearest4")
{
    Image image = Image("../../../../tests/image.png");
    image.height = 3;
    image.width = 3;
    image.pixels = std::vector<glm::vec3> {
        glm::vec3 { 1 }, glm::vec3 { 2 }, glm::vec3 { 3 }, glm::vec3 { 4 }, glm::vec3 { 5 }, glm::vec3 { 6 }, glm::vec3 { 7 }, glm::vec3 { 8 }, glm::vec3 { 9 }
    };
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0, 1 }), glm::vec3 { 1 });
}


TEST_CASE("nearest5")
{
    Image image = Image("../../../../tests/image.png");
    image.height = 3;
    image.width = 3;
    image.pixels = std::vector<glm::vec3> {
        glm::vec3 { 1 }, glm::vec3 { 2 }, glm::vec3 { 3 }, glm::vec3 { 4 }, glm::vec3 { 5 }, glm::vec3 { 6 }, glm::vec3 { 7 }, glm::vec3 { 8 }, glm::vec3 { 9 }
    };
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0.8, 0.5 }), glm::vec3 { 6 });
}

TEST_CASE("nearest6")
{
    Image image = Image("../../../../tests/image.png");
    image.height = 5;
    image.width = 5;
    std::vector<glm::vec3> pixels;

    for (float i = 1; i <= 25; i++) {
        pixels.push_back(glm::vec3{i});
    }
    image.pixels = pixels;

    vec3Check(sampleTextureNearest(image, glm::vec2 { 0.5, 0.5 }), glm::vec3 { 13 });
}


TEST_CASE("nearest7")
{
    Image image = Image("../../../../tests/image.png");
    image.height = 5;
    image.width = 5;
    std::vector<glm::vec3> pixels;

    for (float i = 1; i <= 25; i++) {
        pixels.push_back(glm::vec3 { i });
    }
    image.pixels = pixels;

    vec3Check(sampleTextureNearest(image, glm::vec2 { 0.3, 0.7 }), glm::vec3 { 7 });
}