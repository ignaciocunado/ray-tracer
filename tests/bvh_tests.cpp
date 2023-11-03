// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "shading.h"
#include <limits>

#include "intersect.h"

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()


using Primitive = BVHInterface::Primitive;

Features getFeatures(bool enableBVH = true)
{
    return Features { .enableAccelStructure = enableBVH };
}

Scene getScene(SceneType type = SceneType::CornellBox)
{
    return loadScenePrebuilt(type, DATA_DIR);
}

Scene precalcedScenes[6] = {
    getScene(SceneType::SingleTriangle),
    getScene(SceneType::Cube),
    getScene(SceneType::CornellBox),
    getScene(SceneType::Monkey),
    getScene(SceneType::Teapot),
    getScene(SceneType::Dragon)
};

BVH precalcedBVHs[6] = {
    BVH(precalcedScenes[0], getFeatures()),
    BVH(precalcedScenes[1], getFeatures()),
    BVH(precalcedScenes[2], getFeatures()),
    BVH(precalcedScenes[3], getFeatures()),
    BVH(precalcedScenes[4], getFeatures()),
    BVH(precalcedScenes[5], getFeatures())
};

Scene getPrecalcedScene(SceneType type)
{
    switch (type) {
    case SingleTriangle:
        return precalcedScenes[0];
    case Cube:
        return precalcedScenes[1];
    case CornellBox:
        return precalcedScenes[2];
    case Monkey:
        return precalcedScenes[3];
    case Teapot:
        return precalcedScenes[4];
    case Dragon:
        return precalcedScenes[5];
    default:
        return Scene {}; // dummy Scene just in case
    }
}

BVH getPrecalcedBVH(SceneType type) {
    switch (type) {
    case SingleTriangle:
        return precalcedBVHs[0];
    case Cube:
        return precalcedBVHs[1];
    case CornellBox:
        return precalcedBVHs[2];
    case Monkey:
        return precalcedBVHs[3];
    case Teapot:
        return precalcedBVHs[4];
    case Dragon:
        return precalcedBVHs[5];
    default:
        return BVH(Scene{}, Features{}); // dummy BVH just in case
    }
}

HitInfo createEmptyHitInfo() {
    return HitInfo {
        .normal = glm::vec3(0.0f),
        .barycentricCoord = glm::vec3(0.0f),
        .texCoord = glm::vec2(0.0f),
        .material = Material {
            .kd = glm::vec3(0.0f) }
    };
}

void checkVec2(const glm::vec2& actual, const glm::vec2& expected, const float eps = 1e-9f)
{
    for (int coord = 0; coord < 2; coord++) {
        REQUIRE_THAT(actual[coord], Catch::Matchers::WithinAbs(expected[coord], eps));
    }
}

void checkVec3(const glm::vec3& actual, const glm::vec3& expected, const float eps = 1e-9f) {
    for (int coord = 0; coord < 3; coord++) {
        REQUIRE_THAT(actual[coord], Catch::Matchers::WithinAbs(expected[coord], eps));
    }
}

void checkVertex(const Vertex& actual, const Vertex& expected, const float eps = 1e-9) {
    checkVec3(actual.position, expected.position, eps);
    checkVec3(actual.normal, expected.normal, eps);
    checkVec2(actual.texCoord, expected.texCoord, eps);
}

void checkPrimitive(const Primitive& actual, const Primitive& expected, const float eps = 1e-9f) {
    REQUIRE(actual.meshID == expected.meshID);
    checkVertex(actual.v0, expected.v0, eps);
    checkVertex(actual.v1, expected.v1, eps);
    checkVertex(actual.v2, expected.v2, eps);
}

void checkPrimitiveSpan(const std::span<Primitive> actual, const std::span<Primitive> expected, const float eps = 1e-9f) {
    REQUIRE(actual.size() == expected.size());
    for (int i = 0; i < actual.size(); i++) {
        checkPrimitive(actual[i], expected[i], eps);
    }
}

void checkRays(const Ray& actual, const Ray& expected, const float eps = 1e-9f) {
    checkVec3(actual.origin, expected.origin, eps);
    checkVec3(actual.direction, expected.direction, eps);

    REQUIRE_THAT(actual.t, Catch::Matchers::WithinAbs(expected.t, eps));
}

void checkMaterials(const Material& actual, const Material& expected, const float eps = 1e-9f) {
    checkVec3(actual.kd, expected.kd, eps);
    checkVec3(actual.ks, expected.ks, eps);
    REQUIRE_THAT(actual.shininess,
        Catch::Matchers::WithinAbs(expected.shininess, eps));
    REQUIRE_THAT(actual.transparency,
        Catch::Matchers::WithinAbs(expected.transparency, eps));
    // Might be good idea to check material.kdTexture also. Not needed now.
}

void checkHitInfo(const HitInfo& actual, const HitInfo& expected, const float eps = 1e-9f) {
    checkVec3(actual.barycentricCoord, expected.barycentricCoord, eps);
    checkVec3(actual.normal, expected.normal, eps);
    checkVec2(actual.texCoord, expected.texCoord, eps);
    checkMaterials(actual.material, expected.material, eps);
}

const uint32_t samplerSeed = 123;
Sampler sampler(samplerSeed);

float randomValue(const float from, const float to) {
    return from + (to-from) * sampler.next_1d();
}

glm::vec3 randomVec3() {
    return glm::vec3(
        randomValue(-3.0f, 3.0f),
        randomValue(-3.0f, 3.0f),
        randomValue(-3.0f, 3.0f)
    );
}

Ray randomRay() {
    return Ray {
        .origin = randomVec3(),
        .direction = randomVec3()
    };
}

template <typename T>
std::vector<std::vector<T>> getAllPermutations(const std::vector<T>& original) {
    std::vector<std::vector<T>> result {};

    std::vector<int> permutationOfIndexes(original.size());
    for (int i = 0; i < original.size(); i++) {
        permutationOfIndexes[i] = i;
    }

    do {
        std::vector<T> currentPermutation(original.size());
        for (int i = 0; i < original.size(); i++) {
            currentPermutation[i] = original[permutationOfIndexes[i]];
        }
        result.push_back(currentPermutation);
    } while (std::next_permutation(permutationOfIndexes.begin(), permutationOfIndexes.end()));

    return result;
}

TEST_CASE("Compute AABB for one Primitive")
{
    std::array<Vertex, 3> triangleVertices {
        Vertex { .position = glm::vec3(-2.0f, 2.0f, -3.0f) },
        Vertex { .position = glm::vec3(1.0f, 5.0f, 3.0f) },
        Vertex { .position = glm::vec3(10.0f, -1.0f, 4.0f) }
    };

    AxisAlignedBox expected {
        .lower = glm::vec3(-2.0f, -1.0f, -3.0f),
        .upper = glm::vec3(10.0f, 5.0f, 4.0f)
    };

    std::array<int, 3> permutation { 0, 1, 2 };

    // In this test, all permutations of triangle vertices are checked.
    int permutationId = 0;
    do {
        BVHInterface::Primitive triangle {
            .v0 = triangleVertices[permutation[0]],
            .v1 = triangleVertices[permutation[1]],
            .v2 = triangleVertices[permutation[2]]
        };

        SECTION("computePrimitiveAABB. Permutation #" + std::to_string(permutationId))
        {
            AxisAlignedBox result = computePrimitiveAABB(triangle);
            checkVec3(result.lower, expected.lower);
            checkVec3(result.upper, expected.upper);
        }

        SECTION("computeSpanAABB. Permutation #" + std::to_string(permutationId))
        {
            AxisAlignedBox result = computeSpanAABB(std::vector { triangle });
            checkVec3(result.lower, expected.lower);
            checkVec3(result.upper, expected.upper);
        }

        permutationId++;
    } while (std::next_permutation(permutation.begin(), permutation.end()));
}

TEST_CASE("Compute AABB for Span of Primitives")
{
    std::vector<BVHInterface::Primitive> triangles {
        // triangle 0
        BVHInterface::Primitive {
            .v0 = Vertex { .position = glm::vec3(-2.0f, 2.0f, -3.0f) },
            .v1 = Vertex { .position = glm::vec3(1.0f, 5.0f, 3.0f) },
            .v2 = Vertex { .position = glm::vec3(10.0f, -1.0f, 4.0f) } },

        // triangle 1
        BVHInterface::Primitive {
            .v0 = Vertex { .position = glm::vec3(-50.0f, 4.0f, -8.0f) },
            .v1 = Vertex { .position = glm::vec3(2.0f, 4.0f, 10.0f) },
            .v2 = Vertex { .position = glm::vec3(10.0f, -1.0f, 44.0f) } },

        // triangle 2
        BVHInterface::Primitive {
            .v0 = Vertex { .position = glm::vec3(5.0f, 48.0f, -38.0f) },
            .v1 = Vertex { .position = glm::vec3(1.0f, -17.0f, -3.0f) },
            .v2 = Vertex { .position = glm::vec3(0.0f, -11.0f, 4.0f) } }
    };

    AxisAlignedBox expected {
        .lower = glm::vec3(-50.0f, -17.0f, -38.0f),
        .upper = glm::vec3(10.0f, 48.0f, 44.0f)
    };

    AxisAlignedBox result = computeSpanAABB(triangles);

    checkVec3(result.lower, expected.lower);
    checkVec3(result.upper, expected.upper);
}

TEST_CASE("Compute Centroid of Primitive")
{
    SECTION("Equilateral triangle") {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(0.0f, -0.5f, 0.0f) },
            .v1 = Vertex { .position = glm::vec3(0.0f, 0.5f, 0.0f) },
            .v2 = Vertex { .position = glm::vec3(0.0f, 0.0f, sqrt(3.0f)/2.0f) }
        };

        const glm::vec3 expected { 0.0f, 0.0f, sqrt(3.0f)/6.0f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected);
    }

    SECTION("Degenerate triangle - three points in line") {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(5.0f, 4.0f, 3.0f) },
            .v1 = Vertex { .position = glm::vec3(0.0f, 0.0f, 0.0f) },
            .v2 = Vertex { .position = glm::vec3(-0.5f, -0.4f, -0.3f) }
        };

        const glm::vec3 expected { 1.5f, 1.2f, 0.9f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected, 1e-6f);
    }

    SECTION("Degenerate triangle - three points in the same point") {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(5.0f, 4.0f, 3.0f) },
            .v1 = Vertex { .position = glm::vec3(5.0f, 4.0f, 3.0f) },
            .v2 = Vertex { .position = glm::vec3(5.0f, 4.0f, 3.0f) }
        };

        const glm::vec3 expected { 5.0f, 4.0f, 3.0f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected);
    }

    SECTION("Triangle #0")
    {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(-2.0f, 2.0f, -3.0f) },
            .v1 = Vertex { .position = glm::vec3(1.0f, 5.0f, 3.0f) },
            .v2 = Vertex { .position = glm::vec3(10.0f, -1.0f, 4.0f) }
        };

        const glm::vec3 expected { 3.0f, 2.0f, 4.0f / 3.0f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected);
    }

    SECTION("Triangle #1")
    {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(-50.0f, 4.0f, -8.0f) },
            .v1 = Vertex { .position = glm::vec3(2.0f, 4.0f, 10.0f) },
            .v2 = Vertex { .position = glm::vec3(10.0f, -1.0f, 44.0f) }
        };

        const glm::vec3 expected { -12.0f - 2.0f / 3.0f, 7.0f / 3.0f, 15.0f + 1.0f / 3.0f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected);
    }

    SECTION("Triangle #2")
    {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(5.0f, 48.0f, -38.0f) },
            .v1 = Vertex { .position = glm::vec3(1.0f, -17.0f, -3.0f) },
            .v2 = Vertex { .position = glm::vec3(0.0f, -11.0f, 4.0f) }
        };

        const glm::vec3 expected { 2.0f, 20.0f / 3.0f, -12.0f - 1.0f / 3.0f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected);
    }
}

TEST_CASE("Compute AABB Longest axis")
{
    SECTION("Case #0: All different, longest is x")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(1.0f, 2.0f, 3.0f),
            .upper = glm::vec3(6.0f, 5.0f, 4.0f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 0);
    }

    SECTION("Case #1: All different, longest is y")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(1.0f, 2.0f, 3.0f),
            .upper = glm::vec3(6.0f, 10.0f, 4.0f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 1);
    }

    SECTION("Case #2: All different, longest is z")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(1.0f, 2.0f, 3.0f),
            .upper = glm::vec3(6.0f, 5.0f, 10.0f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 2);
    }

    SECTION("Case #3: x and y are the same and longest")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(1.0f + 1e-6f, 3.0f + 1e-9f, 2.0f),
            .upper = glm::vec3(6.0f + 1e-6f, 8.0f + 1e-9f, 5.0f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 0);
    }

    SECTION("Case #4: x and z are the same and longest")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(1.0f + 1e-6f, 2.0f, 3.0f + 1e-9f),
            .upper = glm::vec3(6.0f + 1e-6f, 5.0f, 8.0f + 1e-9f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 0);
    }

    SECTION("Case #5: y and z are the same and longest")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(2.0f, 1.0f + 1e-6f, 3.0f + 1e-9f),
            .upper = glm::vec3(5.0f, 6.0f + 1e-6f, 8.0f + 1e-9f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 1);
    }

    SECTION("Case #6: All three axes are the same")
    {
        const AxisAlignedBox box {
            .lower = glm::vec3(1.0f + 1e-6f + 1e-9f, 1.0f + 1e-6f, 3.0f + 1e-9f),
            .upper = glm::vec3(6.0f + 1e-9f + 1e-6f, 6.0f + 1e-6f, 8.0f + 1e-9f)
        };
        const uint32_t result = computeAABBLongestAxis(box);
        CHECK(result == 0);
    }
}

TEST_CASE("Split primitives by median") {

    SECTION("One triangle") {
        const size_t expectedSplitPosition = 1;
        std::vector<Primitive> triangles {
            // triangle 0 with centroid at (-5, -3, 6)
            Primitive {
                .meshID = 0,
                .v0 = Vertex { .position = glm::vec3(-10.0f, -3.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(-5.0f, -3.0f, 4.0f) },
                .v2 = Vertex { .position = glm::vec3(0.0f, -3.0f, 7.0f) } }
        };
        const auto allPermutations = getAllPermutations(triangles);
        const AxisAlignedBox aabb = computeSpanAABB(triangles);

        std::vector<Primitive> expectedOrderX { triangles[0] };
        std::vector<Primitive> expectedOrderY { triangles[0] };
        std::vector<Primitive> expectedOrderZ { triangles[0] };

        for (int permutationId = 0; permutationId < allPermutations.size(); permutationId++) {
            SECTION("Permutation #" + std::to_string(permutationId))
            {
                SECTION("Sort by X axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 0, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderX);
                }

                SECTION("Sort by Y axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 1, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderY);
                }

                SECTION("Sort by Z axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 2, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderZ);
                }
            }
        }
    }

    SECTION("Two triangles") { 
        const size_t expectedSplitPosition = 1;
        std::vector<Primitive> triangles {
            // triangle 0 with centroid at (-5, -3, 6)
            Primitive {
                .meshID = 0,
                .v0 = Vertex { .position = glm::vec3(-10.0f, -3.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(-5.0f, -3.0f, 4.0f) },
                .v2 = Vertex { .position = glm::vec3(0.0f, -3.0f, 7.0f) } },

            // triangle 1 with centroid at (2, -20, 8)
            Primitive {
                .meshID = 1,
                .v0 = Vertex { .position = glm::vec3(2.0f, -19.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(2.0f, -20.0f, 7.0f) },
                .v2 = Vertex { .position = glm::vec3(2.0f, -21.0f, 10.0f) } }
        };
        const auto allPermutations = getAllPermutations(triangles);
        const AxisAlignedBox aabb = computeSpanAABB(triangles);

        std::vector<Primitive> expectedOrderX { triangles[0], triangles[1] };
        std::vector<Primitive> expectedOrderY { triangles[1], triangles[0] };
        std::vector<Primitive> expectedOrderZ { triangles[0], triangles[1] };

        for (int permutationId = 0; permutationId < allPermutations.size(); permutationId++) {
            SECTION("Permutation #" + std::to_string(permutationId))
            {
                SECTION("Sort by X axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 0, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderX);
                }

                SECTION("Sort by Y axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 1, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderY);
                }

                SECTION("Sort by Z axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 2, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderZ);
                }
            }
        }

    }

    SECTION("Three triangles") { 
        const size_t expectedSplitPosition = 2;
        std::vector<Primitive> triangles {
            // triangle 0 with centroid at (-5, -3, 6)
            Primitive {
                .meshID = 0,
                .v0 = Vertex { .position = glm::vec3(-10.0f, -3.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(-5.0f, -3.0f, 4.0f) },
                .v2 = Vertex { .position = glm::vec3(0.0f, -3.0f, 7.0f) } },

            // triangle 1 with centroid at (2, -20, 8)
            Primitive {
                .meshID = 1,
                .v0 = Vertex { .position = glm::vec3(2.0f, -19.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(2.0f, -20.0f, 7.0f) },
                .v2 = Vertex { .position = glm::vec3(2.0f, -21.0f, 10.0f) } },

            // triangle 2 with centroid at (10, -1, 5)
            Primitive {
                .meshID = 2,
                .v0 = Vertex { .position = glm::vec3(10.0f, -1.0f, 5.0f) },
                .v1 = Vertex { .position = glm::vec3(10.0f, -1.0f, 5.0f) },
                .v2 = Vertex { .position = glm::vec3(10.0f, -1.0f, 5.0f) } }
        
        };
        const auto allPermutations = getAllPermutations(triangles);
        const AxisAlignedBox aabb = computeSpanAABB(triangles);

        std::vector<Primitive> expectedOrderX { triangles[0], triangles[1], triangles[2] };
        std::vector<Primitive> expectedOrderY { triangles[1], triangles[0], triangles[2] };
        std::vector<Primitive> expectedOrderZ { triangles[2], triangles[0], triangles[1] };

        for (int permutationId = 0; permutationId < allPermutations.size(); permutationId++) {
            SECTION("Permutation #" + std::to_string(permutationId))
            {
                SECTION("Sort by X axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 0, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderX);
                }

                SECTION("Sort by Y axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 1, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderY);
                }

                SECTION("Sort by Z axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 2, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderZ);
                }
            }
        }
    }

    SECTION("Four triangles") {
        const size_t expectedSplitPosition = 2;
        std::vector<Primitive> triangles {
            // triangle 0 with centroid at (-5, -3, 6)
            Primitive {
                .meshID = 0,
                .v0 = Vertex { .position = glm::vec3(-10.0f, -3.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(-5.0f, -3.0f, 4.0f) },
                .v2 = Vertex { .position = glm::vec3(0.0f, -3.0f, 7.0f) } },

            // triangle 1 with centroid at (2, -20, 8)
            Primitive {
                .meshID = 1,
                .v0 = Vertex { .position = glm::vec3(2.0f, -19.0f, 7.0f) },
                .v1 = Vertex { .position = glm::vec3(2.0f, -20.0f, 7.0f) },
                .v2 = Vertex { .position = glm::vec3(2.0f, -21.0f, 10.0f) } },

            // triangle 2 with centroid at (10, -1, 5)
            Primitive {
                .meshID = 2,
                .v0 = Vertex { .position = glm::vec3(10.0f, -1.0f, 5.0f) },
                .v1 = Vertex { .position = glm::vec3(10.0f, -1.0f, 5.0f) },
                .v2 = Vertex { .position = glm::vec3(10.0f, -1.0f, 5.0f) } },

            // triangle 3 with centroid at (0, 12, 7)
            Primitive {
                .meshID = 3,
                .v0 = Vertex { .position = glm::vec3(-1.0f, 13.0f, 5.0f) },
                .v1 = Vertex { .position = glm::vec3(0.0f, 12.0f, 8.0f) },
                .v2 = Vertex { .position = glm::vec3(1.0f, 11.0f, 8.0f) } }
        };
        const auto allPermutations = getAllPermutations(triangles);
        const AxisAlignedBox aabb = computeSpanAABB(triangles);

        std::vector<Primitive> expectedOrderX { triangles[0], triangles[3], triangles[1], triangles[2] };
        std::vector<Primitive> expectedOrderY { triangles[1], triangles[0], triangles[2], triangles[3] };
        std::vector<Primitive> expectedOrderZ { triangles[2], triangles[0], triangles[3], triangles[1] };

        for (int permutationId = 0; permutationId < allPermutations.size(); permutationId++) {
            SECTION("Permutation #" + std::to_string(permutationId))
            {
                SECTION("Sort by X axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 0, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderX);
                }

                SECTION("Sort by Y axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 1, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderY);
                }

                SECTION("Sort by Z axis")
                {
                    std::vector<Primitive> curTriangles = allPermutations[permutationId];
                    size_t resultSplitPosition = splitPrimitivesByMedian(aabb, 2, curTriangles);

                    CHECK(resultSplitPosition == expectedSplitPosition);
                    checkPrimitiveSpan(curTriangles, expectedOrderZ);
                }
            }
        }
    }
}

TEST_CASE("BVH construction - prebuilt scenes")
{
    // Note: for large scenes I cannot be sure that my BVH construction
    // is correct since I cannot write everything by hand on the paper.
    // But at least these tests check that BVH construction does not change
    // when the code is changed. (So, these tests are still quite good to catch new bugs).

    SECTION("Scene: Single Triangle")
    {
        const Scene scene = getPrecalcedScene(SceneType::SingleTriangle);
        const Features features = getFeatures();
        const BVH bvh = getPrecalcedBVH(SceneType::SingleTriangle);

        CHECK(bvh.nodes().size() == 2); // two nodes: Root and Dummy
        CHECK(bvh.primitives().size() == 1);
        CHECK(bvh.numLeaves() == 1);
        CHECK(bvh.numLevels() == 1);
    }

    SECTION("Scene: Cube")
    {
        const Scene scene = getPrecalcedScene(SceneType::Cube);
        const Features features = getFeatures();
        const BVH bvh = getPrecalcedBVH(SceneType::Cube);

        CHECK(bvh.nodes().size() == 8);
        CHECK(bvh.primitives().size() == 12);
        CHECK(bvh.numLeaves() == 4);
        CHECK(bvh.numLevels() == 3);
    }

    SECTION("Scene: Cornell Box")
    {
        const Scene scene = getPrecalcedScene(SceneType::CornellBox);
        const Features features = getFeatures();
        const BVH bvh = getPrecalcedBVH(SceneType::CornellBox);

        CHECK(bvh.nodes().size() == 16);
        CHECK(bvh.primitives().size() == 32);
        CHECK(bvh.numLeaves() == 8);
        CHECK(bvh.numLevels() == 4);
    }

    SECTION("Scene: Monkey")
    {
        const Scene scene = getPrecalcedScene(SceneType::Monkey);
        const Features features = getFeatures();
        const BVH bvh = getPrecalcedBVH(SceneType::Monkey);

        CHECK(bvh.nodes().size() == 512);
        CHECK(bvh.primitives().size() == 967);
        CHECK(bvh.numLeaves() == 256);
        CHECK(bvh.numLevels() == 9);
    }

    SECTION("Scene: Teapot")
    {
        const Scene scene = getPrecalcedScene(SceneType::Teapot);
        const Features features = getFeatures();
        const BVH bvh = getPrecalcedBVH(SceneType::Teapot);

        CHECK(bvh.nodes().size() == 8192);
        CHECK(bvh.primitives().size() == 15704);
        CHECK(bvh.numLeaves() == 4096);
        CHECK(bvh.numLevels() == 13);
    }

    SECTION("Scene: Dragon")
    {
        const Scene scene = getPrecalcedScene(SceneType::Dragon);
        const Features features = getFeatures();
        const BVH bvh = getPrecalcedBVH(SceneType::Dragon);

        CHECK(bvh.nodes().size() == 65536);
        CHECK(bvh.primitives().size() == 87130);
        CHECK(bvh.numLeaves() == 32768);
        CHECK(bvh.numLevels() == 16);
    }
}

TEST_CASE("BVH intersection - prebuilt scenes") {
    // In these tests many random rays are generated and then I test
    // whether the intersection when BVH feature flag is True is the same 
    // as the intersection when BVH flag is False.

    const Features featuresWhenTrue = getFeatures();
    const Features featuresWhenFalse = getFeatures(false);

    const std::vector<SceneType> scenes {
        SceneType::SingleTriangle,
        SceneType::Cube,
        SceneType::CornellBox,
        SceneType::Monkey,
        SceneType::Teapot,
        SceneType::Dragon
    };
    const std::vector<std::string> sceneNames {
        "Single Triangle", "Cube", "Cornell Box", "Monkey", "Teapot", "Dragon"
    };
    const std::vector<int> iterations { 
        5, 20, 100, 300, 300, 300
    };

    for (int sceneId = 0; sceneId < scenes.size(); sceneId++) {
        SECTION("Scene: " + sceneNames[sceneId])
        {
            const BVH bvh = getPrecalcedBVH(scenes[sceneId]);
            const Scene curScene = getPrecalcedScene(scenes[sceneId]);

            RenderState stateWhenTrue {
                .scene = curScene,
                .features = featuresWhenTrue,
                .bvh = bvh
            };

            RenderState stateWhenFalse {
                .scene = curScene, 
                .features = featuresWhenFalse,
                .bvh = bvh
            };

            for (int iteration = 0; iteration < iterations[sceneId]; iteration++) {
                Ray rayWhenTrue = randomRay();
                Ray rayWhenFalse = rayWhenTrue; // copy the same ray
                HitInfo hitInfoWhenTrue = createEmptyHitInfo();
                HitInfo hitInfoWhenFalse = createEmptyHitInfo();

                bool isHitWhenTrue = bvh.intersect(stateWhenTrue, rayWhenTrue, hitInfoWhenTrue);
                bool isHitWhenFalse = bvh.intersect(stateWhenFalse, rayWhenFalse, hitInfoWhenFalse);

                CHECK(isHitWhenTrue == isHitWhenFalse);
                if (isHitWhenTrue) {
                    checkRays(rayWhenTrue, rayWhenFalse);
                    checkHitInfo(hitInfoWhenTrue, hitInfoWhenFalse);
                }
            }

        }
    }
}

TEST_CASE("BVH - pyramid inside another pyramid")
{
    Features featuresWhenTrue = getFeatures();
    Features featuresWhenFalse = getFeatures(false);

    const std::vector<Vertex> vertices {
        Vertex { .position = glm::vec3(1, -2, -2) },
        Vertex { .position = glm::vec3(5, -1, 1) },
        Vertex { .position = glm::vec3(-4, 0, 2) },
        Vertex { .position = glm::vec3(2, 4, 0) },
        Vertex { .position = glm::vec3(0.1, -0.2, -0.2) },
        Vertex { .position = glm::vec3(0.5, -0.1, 0.1) }, 
        Vertex { .position = glm::vec3(-0.4, 0.0, 0.2) }, 
        Vertex { .position = glm::vec3(0.2, 0.4, 0.0) }
    };

    Scene scene {
        .meshes = {
            Mesh { // outside pyramid
                .vertices = { vertices[0], vertices[1], vertices[2], vertices[3] },
                .triangles = { glm::uvec3(0, 1, 2), glm::uvec3(0, 1, 3), glm::uvec3(0, 2, 3), glm::uvec3(1, 2, 3) } },
            Mesh { // inside pyramid (same as outer one, but 10x smaller)
                .vertices = { vertices[4], vertices[5], vertices[6], vertices[7] },
                .triangles = { glm::uvec3(0, 1, 2), glm::uvec3(0, 1, 3), glm::uvec3(0, 2, 3), glm::uvec3(1, 2, 3) } } }
    };

    BVH bvh(scene, featuresWhenTrue);
    RenderState stateWhenTrue {
        .scene = scene,
        .features = featuresWhenTrue,
        .bvh = bvh
    };
    RenderState stateWhenFalse {
        .scene = scene,
        .features = featuresWhenFalse,
        .bvh = bvh
    };

    {
        // Tests for construction
        CHECK(bvh.numLeaves() == 2);
        CHECK(bvh.numLevels() == 2);

        const std::span<BVH::Node> nodes = bvh.nodes();

        // node 0
        checkVec3(nodes[0].aabb.lower, glm::vec3(-4, -2, -2));
        checkVec3(nodes[0].aabb.upper, glm::vec3(5, 4, 2));
        CHECK(!nodes[0].isLeaf());
        CHECK(nodes[0].leftChild() == 2);
        CHECK(nodes[0].rightChild() == 3);

        // node 2
        checkVec3(nodes[2].aabb.lower, glm::vec3(-4, -2, -2));
        checkVec3(nodes[2].aabb.upper, glm::vec3(2, 4, 2));
        CHECK(nodes[2].isLeaf());
        CHECK(nodes[2].primitiveCount() == 4);
        CHECK(nodes[2].primitiveOffset() == 0);

        // node 3
        checkVec3(nodes[3].aabb.lower, glm::vec3(-4, -2, -2));
        checkVec3(nodes[3].aabb.upper, glm::vec3(5, 4, 2));
        CHECK(nodes[3].isLeaf());
        CHECK(nodes[3].primitiveCount() == 4);
        CHECK(nodes[3].primitiveOffset() == 4);

        // check stored primitives
        const std::span<Primitive> primitives = bvh.primitives();
        const std::vector<Primitive> expectedPrimitives {
            Primitive { .meshID = 0, .v0 = vertices[0], .v1=vertices[2], .v2=vertices[3] },
            Primitive { .meshID = 1, .v0 = vertices[4], .v1=vertices[6], .v2=vertices[7] },
            Primitive { .meshID = 1, .v0 = vertices[4], .v1=vertices[5], .v2=vertices[6] },
            Primitive { .meshID = 1, .v0 = vertices[5], .v1=vertices[6], .v2=vertices[7] },
            Primitive { .meshID = 1, .v0 = vertices[4], .v1=vertices[5], .v2=vertices[7] },
            Primitive { .meshID = 0, .v0 = vertices[0], .v1=vertices[1], .v2=vertices[2] },
            Primitive { .meshID = 0, .v0 = vertices[1], .v1=vertices[2], .v2=vertices[3] },
            Primitive { .meshID = 0, .v0 = vertices[0], .v1=vertices[1], .v2=vertices[3] }
        };
        REQUIRE(primitives.size() == 8);
        for (int i = 0; i < 8; i++) {
            checkPrimitive(primitives[i], expectedPrimitives[i]);
        }
    }

    {
        // Intersection test #1: Ray from outside
        Ray ray1 { .origin = glm::vec3(-10, 0, 0), .direction = glm::vec3(1, 0, 0) };
        Ray ray2 = ray1; // copy
        Ray expectedRay = ray1; // copy
        expectedRay.t = 9.199999809265f;

        bool expectedHit = true;

        HitInfo hitInfo1 = createEmptyHitInfo();
        HitInfo hitInfo2 = createEmptyHitInfo();

        bool hit1 = bvh.intersect(stateWhenTrue, ray1, hitInfo1);
        bool hit2 = bvh.intersect(stateWhenFalse, ray2, hitInfo2);

        REQUIRE(expectedHit == hit1);
        REQUIRE(expectedHit == hit2);

        checkRays(ray1, expectedRay);
        checkRays(ray2, expectedRay);
        checkHitInfo(hitInfo1, hitInfo2);
    }

    {
        // Intersection test #2: Ray from between pyramids, should hit the inner one.
        Ray ray1 { .origin = glm::vec3(-1, 0, 0), .direction = glm::vec3(1, 0, 0) };
        Ray ray2 = ray1; // copy
        Ray expectedRay = ray1; // copy
        expectedRay.t = 0.200000092387f;

        bool expectedHit = true;

        HitInfo hitInfo1 = createEmptyHitInfo();
        HitInfo hitInfo2 = createEmptyHitInfo();

        bool hit1 = bvh.intersect(stateWhenTrue, ray1, hitInfo1);
        bool hit2 = bvh.intersect(stateWhenFalse, ray2, hitInfo2);

        REQUIRE(expectedHit == hit1);
        REQUIRE(expectedHit == hit2);

        checkRays(ray1, expectedRay);
        checkRays(ray2, expectedRay);
        checkHitInfo(hitInfo1, hitInfo2);
    }

    {
        // Intersection test #3: Ray from inside inner one.
        Ray ray1 { .origin = glm::vec3(0, 0, 0), .direction = glm::vec3(1, 0, 0) };
        Ray ray2 = ray1; // copy
        Ray expectedRay = ray1; // copy
        expectedRay.t = 0.3249999880790f;

        bool expectedHit = true;

        HitInfo hitInfo1 = createEmptyHitInfo();
        HitInfo hitInfo2 = createEmptyHitInfo();

        bool hit1 = bvh.intersect(stateWhenTrue, ray1, hitInfo1);
        bool hit2 = bvh.intersect(stateWhenFalse, ray2, hitInfo2);

        REQUIRE(expectedHit == hit1);
        REQUIRE(expectedHit == hit2);

        checkRays(ray1, expectedRay);
        checkRays(ray2, expectedRay);
        checkHitInfo(hitInfo1, hitInfo2);
    }

    {
        // Intersection test #4: Ray from between pyramids, should hit the outer one.
        Ray ray1 { .origin = glm::vec3(1, 0, 0), .direction = glm::vec3(1, 0, 0) };
        Ray ray2 = ray1; // copy
        Ray expectedRay = ray1; // copy
        expectedRay.t = 2.25f;

        bool expectedHit = true;

        HitInfo hitInfo1 = createEmptyHitInfo();
        HitInfo hitInfo2 = createEmptyHitInfo();

        bool hit1 = bvh.intersect(stateWhenTrue, ray1, hitInfo1);
        bool hit2 = bvh.intersect(stateWhenFalse, ray2, hitInfo2);

        REQUIRE(expectedHit == hit1);
        REQUIRE(expectedHit == hit2);

        checkRays(ray1, expectedRay);
        checkRays(ray2, expectedRay);
        checkHitInfo(hitInfo1, hitInfo2);
    }

    {
        // Intersection test #5: Ray from outside, no hit.
        Ray ray1 { .origin = glm::vec3(10, 0, 0), .direction = glm::vec3(1, 0, 0) };
        Ray ray2 = ray1; // copy

        bool expectedHit = false;

        HitInfo hitInfo1 = createEmptyHitInfo();
        HitInfo hitInfo2 = createEmptyHitInfo();

        bool hit1 = bvh.intersect(stateWhenTrue, ray1, hitInfo1);
        bool hit2 = bvh.intersect(stateWhenFalse, ray2, hitInfo2);

        REQUIRE(expectedHit == hit1);
        REQUIRE(expectedHit == hit2);
    }
}
