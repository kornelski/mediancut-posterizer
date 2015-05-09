/**
 Median Cut Posterizer
 © 2011-2012 Kornel Lesiński.

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
#include <getopt.h>
#include "png.h"
#include "rwpng.h"

void optimizeForAverageFilter(
    unsigned char pixels[],
    int width, int height,
    int quantization);

#ifndef MAX
 #define MAX(a,b) ((a)>=(b)?(a):(b))
#endif

#define BOTH(a) ((a).color + (a).alpha)

typedef struct {
    unsigned char r,g,b,a;
} rgba_pixel;

typedef struct {
    unsigned int indices[256];
} palette;

typedef struct {
    double color, alpha;
} hist_entry;

inline static void pal_set(palette *pal, const unsigned int val) {
    pal->indices[val] = val;
}

inline static bool pal_isset(const palette *pal, const unsigned int val) {
    return pal->indices[val] == val;
}

static void pal_init(palette *pal) {
    memset(pal->indices, 0, sizeof(pal->indices));
}

static void interpolate_palette_front(const palette *pal, unsigned int mapping[], const bool dither);
static void voronoi(const hist_entry histogram[static 256], palette *pal);
static double palette_error(const hist_entry histogram[static 256], const palette *palette_orig);
static void interpolate_palette_back(const palette *pal, unsigned int mapping[]);
static void posterize(png24_image *img, unsigned int maxlevels, const double maxerror, bool dither, bool verbose);

inline static double int_to_linear(unsigned int value)
{
    return value/255.0;
}

// *256 is not off-by-one error.
inline static unsigned int linear_to_int(const double value)
{
    const double g = value*256.0;
    return g < 255.0 ? g : 255;
}

static double image_gamma, gamma_lut[256];
static void set_gamma(const double gamma)
{
    image_gamma = gamma;
    for(int i=0; i < 256; i++) gamma_lut[i] = pow(int_to_linear(i), gamma);
}

// Converts gamma 2.2 to linear unit value. Linear color is required for preserving brightness (esp. when dithering).
// call set_gamma first
inline static double gamma_to_linear(unsigned int value)
{
    return gamma_lut[value];
}

// Reverses gamma_to_linear.
inline static unsigned int linear_to_gamma(const double value)
{
    return linear_to_int(pow(value, 1.0/image_gamma));
}

// median cut "box" in this implementation is actually a line,
// since it only needs to track lowest/highest intensity
struct box {
    double sum, variance;
    unsigned int start, end;
};

// helper function that gives integer intensity (palette index) from given weights.
// NB: in this function color is linear 0..1, alpha is 0..255!
inline static unsigned int index_from_weights(hist_entry weight, hist_entry sum)
{
    const double color_gamma = weight.color ? linear_to_gamma(sum.color/weight.color) * weight.color : 0;
    const double mixed_linear = (color_gamma + sum.alpha) / (BOTH(weight) * 255.0);
    return linear_to_int(mixed_linear);
}

// average values in a "box" proportionally to frequency of their occurence
// returns linear value (which is a mix of color and alpha components, so can't be gamma-corrected later)
static double weighted_avg_linear(const unsigned int start, const unsigned int end, const hist_entry histogram[static 256])
{
    double weight=0,sum=0;
    for(unsigned int val=start; val < end; val++) {
        weight += BOTH(histogram[val]);
        sum += gamma_to_linear(val)*histogram[val].color + int_to_linear(val)*histogram[val].alpha;
    }
    return weight ? sum/weight : 0;
}

// returns integer index that from weighed average and applies gamma correction proportionally to amount of color
static unsigned int weighted_avg_int(const unsigned int start, const unsigned int end, const hist_entry histogram[static 256])
{
    hist_entry weight = {0};
    hist_entry sum = {0};

    for(unsigned int val=start; val < end; val++) {
        weight.color += histogram[val].color;
        weight.alpha += histogram[val].alpha;
        sum.color += histogram[val].color * gamma_to_linear(val);
        sum.alpha += histogram[val].alpha * val;
    }

    return index_from_weights(weight, sum);
}

// variance (AKA second moment) of the box. Measures how much "spread" the values are
static double variance_in_range(const unsigned int start, const unsigned int end, const hist_entry histogram[static 256])
{
    const double avg = weighted_avg_linear(start, end, histogram);

    double sum=0;
    for(unsigned int val=start; val < end; val++) {
        const double color_delta = avg-gamma_to_linear(val);
        const double alpha_delta = avg-int_to_linear(val);
        sum += color_delta*color_delta*histogram[val].color;
        sum += alpha_delta*alpha_delta*histogram[val].alpha;
    }
    return sum;
}

static double variance(const struct box box, const hist_entry histogram[static 256])
{
    return variance_in_range(box.start, box.end, histogram);
}

// Square error. Estimates how well palette "fits" the histogram.
static double palette_error(const hist_entry histogram[static 256], const palette *pal)
{
    unsigned int mapping[256];

    // the input palette has gaps
    interpolate_palette_front(pal, mapping, false);

    double sum=0, px=0;
    for (unsigned int i=0; i < 256; i++) {
        double color_delta = gamma_to_linear(i)-gamma_to_linear(mapping[i]);
        double alpha_delta = int_to_linear(i)-int_to_linear(mapping[i]);
        sum += color_delta*color_delta*histogram[i].color;
        sum += alpha_delta*alpha_delta*histogram[i].alpha;
        px += BOTH(histogram[i]);
    }
    return sum/px;
}

// converts boxes to palette.
// palette here is a sparse array where elem[x]=x is taken, elem[x]=0 is free (except x=0)
static void palette_from_boxes(const struct box boxes[], const int numboxes, const hist_entry histogram[static 256], palette *pal)
{
    pal_init(pal);

    for(int box=0; box < numboxes; box++) {
        pal_set(pal, weighted_avg_int(boxes[box].start, boxes[box].end, histogram));
    }
    pal_set(pal, 0);
    pal_set(pal, 255);
}

/*
 1-dimensional median cut, using variance for "largest" box
*/
static unsigned int reduce(const unsigned int maxlevels, const double maxerror, const hist_entry histogram[static 256], palette *pal)
{
    unsigned int numboxes=1;
    struct box boxes[256];

    // build the first "box" that encompasses all values
    boxes[0].start=1; // skip first and last entry, as they're always included
    boxes[0].end=255;
    boxes[0].sum=0;
    for(unsigned int i=boxes[0].start; i < boxes[0].end; i++) boxes[0].sum += BOTH(histogram[i]);
    boxes[0].variance = 1; // irrelevant for first box

    while(numboxes < maxlevels) {
        int boxtosplit=-1;
        double largest=0;
        // pick box to split by choosing one with highest variance
        for(int box=0; box < numboxes; box++) {
            if (boxes[box].variance > largest && (boxes[box].end-boxes[box].start)>=2) {
                largest = boxes[box].variance;
                boxtosplit=box;
            }
        }
        if (boxtosplit < 0) {
            break;
        }

        // divide equally by variance
        unsigned int bestsplit=0;
        double minvariance = INFINITY;
        for(unsigned int val=boxes[boxtosplit].start+1; val < boxes[boxtosplit].end-1; val++) {
            const double variance = variance_in_range(boxes[boxtosplit].start, val, histogram)
                                  + variance_in_range(val, boxes[boxtosplit].end, histogram);
            if (variance < minvariance) {
                minvariance = variance;
                bestsplit = val;
            }
        }

        double sum=0;
        for(unsigned int i=boxes[boxtosplit].start; i < bestsplit; i++) sum += BOTH(histogram[i]);

        // create new boxes from halves
        boxes[numboxes].start = boxes[boxtosplit].start;
        boxes[numboxes].end = bestsplit;
        boxes[numboxes].sum = sum;
        boxes[numboxes].variance = variance(boxes[numboxes], histogram);
        boxes[boxtosplit].start = bestsplit;
        boxes[boxtosplit].sum -= boxes[numboxes].sum;
        boxes[boxtosplit].variance = variance(boxes[boxtosplit], histogram);
        numboxes++;

        if (maxerror > 0 && maxerror != INFINITY) {
            palette_from_boxes(boxes, numboxes, histogram, pal);

            voronoi(histogram, pal);

            if (palette_error(histogram, pal) < maxerror) {
                return numboxes;
            }
        }
    }

    palette_from_boxes(boxes, numboxes, histogram, pal);

    return numboxes;
}

// palette1/2 is for even/odd pixels, allowing very simple "ordered" dithering
static void remap(png24_image *img, const palette *pal, bool dither)
{
    unsigned int mapping1[256], mapping2[256];

    if (dither) {
        // front to back. When dithering, it's biased towards nextval
        interpolate_palette_front(pal, mapping1, true);

        // back to front, so dithering bias is the other way.
        interpolate_palette_back(pal, mapping2);
    } else {
        interpolate_palette_front(pal, mapping1, false);
        memcpy(mapping2, mapping1, sizeof(mapping2));
    }

    for(unsigned int i=0; i < img->height; i++) {
        rgba_pixel *const row = (rgba_pixel*)img->row_pointers[i];
        for(unsigned int j=0; j < img->width; j++) {
            const unsigned int *map = (i^j)&1 ? mapping1 : mapping2;
            const rgba_pixel px = row[j];
            if (map[px.a]) {
                row[j] = (rgba_pixel){
                  .r = map[px.r],
                  .g = map[px.g],
                  .b = map[px.b],
                  .a = map[px.a],
                };
            } else {
                // clear "dirty alpha"
                row[j] = (rgba_pixel){0,0,0,0};
            }
        }
    }
}

// it doesn't count unique colors, only intensity values of all channels
static void intensity_histogram(const png24_image *img, hist_entry histogram[static 256])
{
    for(unsigned int i=0; i < img->height; i++) {
        const rgba_pixel *const row = (rgba_pixel*)img->row_pointers[i];
        for(unsigned int j=0; j < img->width; j++) {
            const rgba_pixel px = row[j];
            // opaque colors get more weight
            const double weight = px.a/255.0;

            // color and alpha are tracked separately, because
            // difference between colors is non-linear (gamma applies)
            // e.g. dark colors are less visually distinct than low alpha values
            histogram[px.r].color += weight*0.975;
            histogram[px.g].color += weight*0.975;
            histogram[px.b].color += weight*0.975;
            // a little weight in non-gamma values, this is a fudge to fix some banding errors
            histogram[px.r].alpha += weight*0.025;
            histogram[px.g].alpha += weight*0.025;
            histogram[px.b].alpha += weight*0.025;
            histogram[px.a].alpha += 1.0 + 3.0*(1.0-weight);
        }
    }
}

// interpolates front-to-back. If dither is true, it will bias towards one side
static void interpolate_palette_front(const palette *pal, unsigned int mapping[], const bool dither)
{
    unsigned int nextval=0, lastval=0;
    assert(pal_isset(pal,0));
    assert(pal_isset(pal,255));

    for(unsigned int val=0; val < 256; val++) {
        if (pal_isset(pal, val)) {
            lastval = val;
            for(unsigned int j=val+1; j < 256; j++) {
                if (pal_isset(pal, j)) {nextval=j; break;}
            }
        }
        const double lastvaldiff = (int_to_linear(val) - int_to_linear(lastval));
        const double nextvaldiff = (int_to_linear(nextval) - int_to_linear(val));
        if (!dither) {
            mapping[val] = lastvaldiff < nextvaldiff ? lastval : nextval;
        } else {
            mapping[val] = lastvaldiff/2 < nextvaldiff ? lastval : nextval;
        }
    }
}

// interpolates back-to-front. Always biased for dither.
static void interpolate_palette_back(const palette *pal, unsigned int mapping[])
{
    unsigned int nextval=255, lastval=255;

    for(int val=255; val >= 0; val--) {
        if (pal_isset(pal, val)) {
            lastval = val;
            for(int j=val-1; j >= 0; j--) {
                if (pal_isset(pal, j)) {nextval=j; break;}
            }
        }
        const double lastvaldiff = (int_to_linear(val) - int_to_linear(lastval));
        const double nextvaldiff = (int_to_linear(nextval) - int_to_linear(val));
        mapping[val] = lastvaldiff/2 >= nextvaldiff ? lastval : nextval;
    }
}

static void usage(const char *exepath)
{
    const char *name = strrchr(exepath, '/');
    if (name) name++; else name = exepath;
    fprintf(stderr, "Median Cut PNG Posterizer 2.1 (2015).\n" \
    "Usage: %s [-vdb] [-Q <quality>] [levels] [input file] [output file]\n\n" \
    "Specify number of levels (2-255) or quality (10-100).\n" \
    "-b blurize mode (uses diagonal averaging filter, recommended)\n" \
    "-d enables dithering\n" \
    "-v verbose output (to stderr)\n\n" \
    "If files are not specified stdin and stdout is used.\n"
    "%s -Q 95 in.png out.png\n", name, name);
}

// performs voronoi iteration (mapping histogram to palette and creating new palette from remapped values)
// this shifts palette towards local optimum
static void voronoi(const hist_entry histogram[static 256], palette *pal)
{
    unsigned int mapping[256];

    interpolate_palette_front(pal, mapping, false);

    hist_entry weights[256] = {{0}};
    hist_entry sums[256] = {{0}};

    // remap palette
    for (unsigned int val=0; val < 256; val++) {
        int best = mapping[val];
        if (0==best || 255==best) continue; // those two are guaranteed to be present, so ignore their influence
        weights[best].color += histogram[val].color;
        weights[best].alpha += histogram[val].alpha;
        sums[best].color += histogram[val].color * gamma_to_linear(val);
        sums[best].alpha += histogram[val].alpha * val;
    }

    pal_init(pal);

    // rebuild palette from remapped averages
    for(unsigned int i=1; i < 255; i++) {
        if (BOTH(weights[i])) {
            pal_set(pal, index_from_weights(weights[i], sums[i]));
        }
    }
    pal_set(pal, 0);
    pal_set(pal, 255);
}


static double quality_to_mse(long quality)
{
    if (quality == 0) return INFINITY;

    // curve fudged to be roughly similar to quality of libjpeg
    // except lowest 10 for really low number of colors
    const double extra_low_quality_fudge = MAX(0,0.016/(0.001+quality) - 0.001);
    return (extra_low_quality_fudge + 2.5/pow(210.0 + quality, 1.2) * (100.1-quality)/100.0) / 6.0;
}

static unsigned int mse_to_quality(double mse)
{
    for(int i=100; i > 0; i--) {
        if (mse <= quality_to_mse(i)) return i;
    }
    return 0;
}

#include <unistd.h>

#if defined(WIN32) || defined(__WIN32__)
#include <fcntl.h>
#include <io.h>
#else
#define setmode(what,ever)
#endif

int main(int argc, char *argv[])
{
    bool dither = false, verbose = false;
    bool blurize = false;
    int quality = 0;

    int ch;
    while ((ch = getopt(argc, argv, "hvdq:Q:b")) != -1) {
        switch (ch) {
            case 'b': blurize = true; break;
            case 'd': dither = true; break;
            case 'v': verbose = true; break;
            case 'q':
            case 'Q':
                quality = atol(optarg);
                break;
            case '?': case 'h':
            default:
                usage(argv[0]);
                return 1;
        }
    }
    int argn = optind;

    int maxlevels = quality > 0 ? 255 : 0;

    if (argn < argc) {
        char *levels_end;
        unsigned long levels = strtoul(argv[argn], &levels_end, 10);
        if (levels_end != argv[argn] && '\0' == levels_end[0]) {
            maxlevels = levels;
            argn++;
        }
    }

    if (maxlevels < 2 || maxlevels > 255) {
        usage(argv[0]);
        return 1;
    }

    FILE *input = stdin;
    const char *input_name = "stdin";
    if (argn < argc && 0 != strcmp("-",argv[argn])) {
        input_name = argv[argn++];
        input = fopen(input_name, "rb");
    }

    FILE *output = stdout;
    const char *output_name = "stdout";
    if (argn < argc && 0 != strcmp("-",argv[argn])) {
        output_name = argv[argn++];
        output = fopen(output_name, "wb");
    }

    if (argn != argc) {
        usage(argv[0]);
        return 1;
    }

    setmode(1, O_BINARY);
    setmode(0, O_BINARY);

    png24_image img;
    pngquant_error retval;

    if ((retval = rwpng_read_image24(input, &img, verbose))) {
        fprintf(stderr, "Error: cannot read PNG from %s\n", input_name);
        return retval;
    }
    if (input != stdin) fclose(input);

    set_gamma(1.0/img.gamma);

    if (blurize) {
        const int quantization = quality ? (103 - quality)/2.6 : 256 - maxlevels;
        optimizeForAverageFilter(img.rgba_data, img.width, img.height, quantization);
    } else {
        double maxerror = quality_to_mse(quality);
        posterize(&img, maxlevels, maxerror, dither, verbose);
    }

    if ((retval = rwpng_write_image24(output, &img, blurize ? PNG_FILTER_VALUE_AVG : PNG_FILTER_VALUE_NONE))) {
        fprintf(stderr, "Error: cannot write PNG to %s\n", output_name);
        return retval;
    }
    if (output != stdout) fclose(output);

    return 0;
}


static void posterize(png24_image *img, unsigned int maxlevels, const double maxerror, bool dither, bool verbose)
{
    hist_entry histogram[256]={{0}};
    intensity_histogram(img, histogram);

    // reserve colors for black and white
    // and omit them from histogram to avoid confusing median cut
    unsigned int reservedcolors=0;
    if (BOTH(histogram[0]) >= 1.0 && maxlevels > 2) {
        maxlevels--;reservedcolors++;
        histogram[0]=(hist_entry){0,0};
    }
    if (BOTH(histogram[255]) >= 1.0 && maxlevels > 2) {
        maxlevels--;reservedcolors++;
        histogram[255]=(hist_entry){0,0};
    }

    palette pal;
    unsigned int levels = reduce(maxlevels, maxerror, histogram, &pal);

    double last_err = INFINITY;
    for(unsigned int j=0; j < 100; j++) {
        voronoi(histogram, &pal);

        double new_err = palette_error(histogram, &pal);
        if (new_err == last_err) break;
        last_err = new_err;
    }

    if (verbose) {
        fprintf(stderr, "MSE=%.3f (Q=%d, %u levels)\n", last_err*65536.0, mse_to_quality(last_err), levels+reservedcolors);
    }

    remap(img, &pal, dither);
}
