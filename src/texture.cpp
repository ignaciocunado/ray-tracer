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

    //Calculate the corresponding index and round down
    float index = floor((j * image.width + i));
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

    // Calculate i and j values
    float i = texCoord[0] * image.width;
    float j = (1 - texCoord[1]) * image.height;
    
    // Interpolate in the x direction
    float i1 = floor(i);
    float i2 = ceil(i);
    if (i2 > image.width) {
        i1--;
        i2--;
    }

    float j1 = floor(j);
    float j2 = ceil(j);
    if (j2 > image.height) {
        j1--;
        j2--;
    }

    float distI1 = abs(i - i1 + 0.5f);
    float distI2 = abs(i2 - i + 0.5f);

    // Interpolate the first row in the x direction
    glm::vec3 pix1 = image.pixels[floor((j1 * image.width + i1))];
    glm::vec3 pix2 = image.pixels[floor((j1 * image.width + i2))];
    glm::vec3 interpolatedX1 = distI2 * pix1 + distI1 * pix2;

    // Interpolate the second row in the x direction
    glm::vec3 pix3 = image.pixels[floor((j2 * image.width + i1))];
    glm::vec3 pix4 = image.pixels[floor((j2 * image.width + i2))];
    glm::vec3 interpolatedX2 = distI2 * pix3 + distI1 * pix4;

    // Interpolate in the y direction
    float distJ1 = abs(j - j1 + 0.5f);
    float distJ2 = abs(j2 - j + 0.5f);

    glm::vec3 res = distJ2 * interpolatedX1 + distJ1 * interpolatedX2;
    
    return res;
}