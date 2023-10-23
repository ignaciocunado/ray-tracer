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

Image imageGenerator5x5()
{
    Image image = Image("../../../../tests/image.png");
    image.height = 5;
    image.width = 5;
    std::vector<glm::vec3> pixels;

    for (float i = 1; i <= 25; i++) {
        pixels.push_back(glm::vec3 { i });
    }
    image.pixels = pixels;

    return image;
}

Image imageGenerator3x3()
{
    Image image = Image("../../../../tests/image.png");
    image.height = 3;
    image.width = 3;
    image.pixels = std::vector<glm::vec3> {
        glm::vec3 { 1 }, glm::vec3 { 2 }, glm::vec3 { 3 }, glm::vec3 { 4 }, glm::vec3 { 5 }, glm::vec3 { 6 }, glm::vec3 { 7 }, glm::vec3 { 8 }, glm::vec3 { 9 }
    };

    return image;
}

TEST_CASE("nearest1")
{
    Image image = imageGenerator3x3();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0 }), glm::vec3 { 7 });
}

TEST_CASE("nearest2")
{
    Image image = imageGenerator3x3();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 1, 0 }), glm::vec3 { 9 });
}

TEST_CASE("nearest3")
{
    Image image = imageGenerator3x3();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 1, 1 }), glm::vec3 { 3 });
}

TEST_CASE("nearest4")
{
    Image image = imageGenerator3x3();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0, 1 }), glm::vec3 { 1 });
}

TEST_CASE("nearest5")
{
    Image image = imageGenerator3x3();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0.8f, 0.5f }), glm::vec3 { 6 });
}

TEST_CASE("nearest6")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0.5f, 0.5f }), glm::vec3 { 13 });
}

TEST_CASE("nearest7")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureNearest(image, glm::vec2 { 0.3f, 0.7f }), glm::vec3 { 7 });
}

TEST_CASE("bilinear1")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0, 0 }), glm::vec3 { 21 });
}

TEST_CASE("bilinear2")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0, 1 }), glm::vec3 { 1 });
}

TEST_CASE("bilinear3")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 1, 0 }), glm::vec3 { 25 });
}

TEST_CASE("bilinear4")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 1, 1 }), glm::vec3 { 5 });
}

TEST_CASE("bilinear5")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.2f, 0.2f }), glm::vec3 { 19 });
}

TEST_CASE("bilinear6")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.1f, 0.1f }), glm::vec3 { 21 });
}

TEST_CASE("bilinear7")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.5f, 0.0f }), glm::vec3 { 23 });
}

TEST_CASE("bilinear8")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.25f, 0.15f }), glm::vec3 { 20.5f });
}

TEST_CASE("bilinear9")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.65f, 0.45f }), glm::vec3 { 15 });
}

TEST_CASE("bilinear10")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.9f, 0.5f }), glm::vec3 { 15 });
}

TEST_CASE("bilinear11")
{
    Image image = imageGenerator5x5();
    vec3Check(sampleTextureBilinear(image, glm::vec2 { 0.01f, 0.5f }), glm::vec3 { 11 });
}