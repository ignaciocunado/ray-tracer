// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "shading.h"
#include <limits>

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()


void checkVec3(const glm::vec3 actual, const glm::vec3 expected, const float eps = 1e-9) {
    for (int coord = 0; coord < 3; coord++) {
        CHECK_THAT(actual[coord], Catch::Matchers::WithinAbs(expected[coord], eps));
    }
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

        checkVec3(result, expected, 1e-6);
    }

    SECTION("Degenerate triangle - three points in line") {
        BVHInterface::Primitive triangle {
            .v0 = Vertex { .position = glm::vec3(5.0f, 4.0f, 3.0f) },
            .v1 = Vertex { .position = glm::vec3(0.0f, 0.0f, 0.0f) },
            .v2 = Vertex { .position = glm::vec3(-0.5f, -0.4f, -0.3f) }
        };

        const glm::vec3 expected { 1.5f, 1.2f, 0.9f };
        const glm::vec3 result = computePrimitiveCentroid(triangle);

        checkVec3(result, expected, 1e-6);
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
