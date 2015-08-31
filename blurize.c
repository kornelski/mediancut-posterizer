/**
 Blurizer
 © 2013 Kornel Lesiński.
 Based on algorithms by Michael Vinther and William MacKay.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details:
 <http://www.gnu.org/copyleft/gpl.html>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include "png.h"
#include "rwpng.h"

typedef struct {
    unsigned char r,g,b,a;
} rgba_pixel;

typedef int colorDelta[4];

static void diffuseColorDeltas(colorDelta **colorError, int x, int *delta);

void optimizeForAverageFilter(
    unsigned char pixels[],
    int width, int height,
    int quantization) {

    int halfStep = quantization / 2;


    const int filterWidth = 5;
    const int filterCenter = 2;
    const int bytesPerPixel = 4;
    const int errorRowCount = 3;
    const int ditheringPrecision = 256;
    int stride = width * bytesPerPixel;
    colorDelta row0[width + filterWidth - 1], row1[width + filterWidth - 1], row2[width + filterWidth - 1];
    colorDelta *colorError[3] = { row0, row1, row2 }, *colorTemp;

    memset(row0, 0, sizeof(row0));
    memset(row1, 0, sizeof(row1));
    memset(row2, 0, sizeof(row2));

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
            colorDelta diffusion = {0,0,0,0};
            diffuseColorDeltas(colorError, x + filterCenter, diffusion);
            for(int c = 0; c < bytesPerPixel; c++) {
                int offset = y*stride + x*bytesPerPixel + c;
                int here = pixels[offset];
                int errorHere = 0;
                if (here > 0 && here < 255) {
                    int up=0, left=0;
                    if (y > 0) {
                        up = pixels[offset-stride];
                    }
                    if (x > 0) {
                        left = pixels[offset-bytesPerPixel];
                    }
                    int average = (up + left) / 2; // PNG average filter

                    int newValue = diffusion[c]/ditheringPrecision + here - average;
                    newValue += halfStep;
                    newValue -= newValue % quantization;
                    newValue += average;
                    if (newValue >= 0 && newValue <= 255) {
                        pixels[offset] = newValue;
                        errorHere = here - newValue;
                    }
                    colorError[0][x + filterCenter][c] = errorHere * ditheringPrecision;
                }
            }
        }
        colorTemp = colorError[errorRowCount-1];
        for(int i = errorRowCount-1; i > 0; i--) {
            colorError[i] = colorError[i-1];
        }
        colorError[0] = colorTemp;
    }
}

static void diffuseColorDeltas(colorDelta **colorError, int x, int *delta) {
    // Sierra dithering
    for( int i = 0; i < 4; i++ ){
        delta[i] += 2 * colorError[2][x-1][i];
        delta[i] += 3 * colorError[2][x][i];
        delta[i] += 2 * colorError[2][x+1][i];
        delta[i] += 2 * colorError[1][x-2][i];
        delta[i] += 4 * colorError[1][x-1][i];
        delta[i] += 5 * colorError[1][x][i];
        delta[i] += 4 * colorError[1][x+1][i];
        delta[i] += 2 * colorError[1][x+2][i];
        delta[i] += 3 * colorError[0][x-2][i];
        delta[i] += 5 * colorError[0][x-1][i];
        if (delta[i] < 0) {
            delta[i] -= 16;
        } else {
            delta[i] += 16;
        }
        delta[i] /= 32;
    }
}

