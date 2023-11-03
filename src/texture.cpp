#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image

    // Calculate i and j values 
    float i = texCoord[0] * image.width;
    float j = (1- texCoord[1]) * image.height;

    if (i == image.width) {
        i--;
    }
    if (j == image.height) {
        j--;
    }

    //Calculate the corresponding index and round down
    float downJ = floor(j);
    float downI = floor(i);
    float index = downJ * image.width + downI;
    return image.pixels[index];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image

    // Calculate coordinates

    float x = texCoord[0] * image.width;
    float y = (1 - texCoord[1]) * image.height;

    // Check if we were given any of the corners or edges
    if (x == 0 && y == 0 || x == image.width && y == 0 || x == 0 && y == image.height || x == image.width && y == image.height || x <= 0.5f || x >= image.width - 0.5f || y <= 0.5f || y >= image.height - 0.5f) {
        if (x >= image.width) {
            x--;
        }
        if (y >= image.height) {
            y--;
        }
        float index = floor(y) * image.width + floor(x);
        return image.pixels[index]; 
    }
    
    // Interpolate in the x direction
    float i1 = round(x) - 0.5f;
    float i2 = round(x) + 0.5f;

    float j1 = round(y) - 0.5f;
    float j2 = round(y) + 0.5f;

    float distI1 = abs(x - i1);
    float distI2 = abs(i2 - x);

    // Interpolate the first row in the x direction
    glm::vec3 pix1 = image.pixels[floor(j1) * image.width + floor(i1)];
    glm::vec3 pix2 = image.pixels[floor(j1) * image.width + floor(i2)];
    glm::vec3 interpolatedX1 = distI2 * pix1 + distI1 * pix2;

    // Interpolate the second row in the x direction
    glm::vec3 pix3 = image.pixels[floor(j2) * image.width + floor(i1)];
    glm::vec3 pix4 = image.pixels[floor(j2) * image.width + floor(i2)];
    glm::vec3 interpolatedX2 = distI2 * pix3 + distI1 * pix4;

    // Interpolate in the y direction
    float distJ1 = abs(y - j1);
    float distJ2 = abs(j2 - y);

    glm::vec3 res = distJ2 * interpolatedX1 + distJ1 * interpolatedX2;
    
    return res;
}