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
#include "png.h"
#include "rwpng.h"

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
static void posterize(png24_image *img, unsigned int maxcolors, const double maxerror, bool dither, bool verbose);

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

// Converts gamma 2.2 to linear unit value. Linear color is required for preserving brightness (esp. when dithering).
double image_gamma = 2.2;
inline static double gamma_to_linear(unsigned int value)
{
    return pow(int_to_linear(value), image_gamma);
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
static unsigned int reduce(const unsigned int maxcolors, const double maxerror, const hist_entry histogram[static 256], palette *pal)
{
    unsigned int numboxes=1;
    struct box boxes[256];

    // build the first "box" that encompasses all values
    boxes[0].start=1; // skip first and last entry, as they're always included
    boxes[0].end=255;
    boxes[0].sum=0;
    for(unsigned int i=boxes[0].start; i < boxes[0].end; i++) boxes[0].sum += BOTH(histogram[i]);
    boxes[0].variance = 1; // irrelevant for first box

    while(numboxes < maxcolors) {
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

        if (maxerror > 0) {
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
            histogram[px.r].color += weight;
            histogram[px.g].color += weight;
            histogram[px.b].color += weight;
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
    fprintf(stderr, "Median Cut PNG Posterizer 1.6 (2013).\n" \
    "Usage: %s [-vd] [-Q <quality>] [levels]\n\n" \
    "Specify number of levels (2-255) or quality (10-100).\n" \
    "-d enables dithering\n" \
    "-v verbose output (to stderr)\n\n" \
    "Image is always read from stdin and written to stdout.\n"
    "%s -d 16 < in.png > out.png\n", name, name);
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

static unsigned int colordifference_ch(const int x, const int y, const int alphas)
{
    // maximum of channel blended on white, and blended on black
    // premultiplied alpha and backgrounds 0/1 shorten the formula
    const int black = x-y, white = black+alphas;
    return black*black + white*white;
}

static unsigned int colordifference(const rgba_pixel px, const rgba_pixel py)
{
    const int alphas = py.a-px.a;
    return colordifference_ch((px.r * px.a)>>8, (py.r * py.a)>>8, alphas) +
           colordifference_ch((px.g * px.a)>>8, (py.g * py.a)>>8, alphas) +
           colordifference_ch((px.b * px.a)>>8, (py.b * py.a)>>8, alphas);
}

static void get_pixel_differences(const rgba_pixel *const  row, const int width, unsigned int *diffs)
{
    // 5 pixels around center, pixels outside the edge replicate 0th pixel
    const rgba_pixel x0 = row[0];
    const rgba_pixel x3 = row[1];
    const rgba_pixel x4 = row[2];
    unsigned int acc[4] = {
        x0.r*3 + x3.r + x4.r,
        x0.g*3 + x3.g + x4.g,
        x0.b*3 + x3.b + x4.b,
        x0.a*3 + x3.a + x4.a,
    };

    for(int x=0; x < width; x++) {
        diffs[x] = colordifference(row[x], (rgba_pixel){
            .r = acc[0] / 5,
            .g = acc[1] / 5,
            .b = acc[2] / 5,
            .a = acc[3] / 5,
        });

        const rgba_pixel prev = row[x > 2 ? x-2 : 0];
        acc[0] -= prev.r;
        acc[1] -= prev.g;
        acc[2] -= prev.b;
        acc[3] -= prev.a;

        const rgba_pixel next = row[x < width-3 ? x+2 : width-1];
        acc[0] += next.r;
        acc[1] += next.g;
        acc[2] += next.b;
        acc[3] += next.a;
    }
}

static void blur_pixel_differences(int width, unsigned int *diffs)
{
    unsigned int prev = diffs[0];
    for(int x=0; x < width; x++) {
        const unsigned int tmp = diffs[x];
        const unsigned int next = diffs[x < width-2 ? x+1 : width-1];
        diffs[x] += prev + next; // no need to divide
        prev = tmp;
    }
}

static void rle_line(rgba_pixel *const row, unsigned int width, unsigned int *diffs, const unsigned int maxerror, const unsigned int mapping1[], const unsigned int mapping2[]);

static void rle(png24_image *img, const double unit_maxerror, const palette *pal, bool dither)
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


    unsigned int maxerror = unit_maxerror * 6.0 * 65536.0; // scale to colordifference output
    unsigned int diffs[img->width];

    for(int y=0; y < img->height; y++) {
        rgba_pixel *const row = (rgba_pixel*)img->row_pointers[y];

        get_pixel_differences(row, img->width, diffs);
        blur_pixel_differences(img->width, diffs);
        blur_pixel_differences(img->width, diffs);
        blur_pixel_differences(img->width, diffs);
        blur_pixel_differences(img->width, diffs);
        blur_pixel_differences(img->width, diffs);


        rle_line((rgba_pixel*)img->row_pointers[y], img->width, diffs, maxerror, mapping1, mapping2);
    }
}

static void remap_line(rgba_pixel *const row, unsigned int width, const unsigned int mapping1[], const unsigned int mapping2[])
{
    for(unsigned int j=0; j < width; j++) {
        const unsigned int *map = ((unsigned int)row ^ j)&1 ? mapping1 : mapping2;// FIXME needs y-offset :(
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

static void rle_line(rgba_pixel *const row, unsigned int width, unsigned int *diffs, const unsigned int maxerror, const unsigned int mapping1[], const unsigned int mapping2[])
{
    if (width < 16) {
        remap_line(row, width, mapping1, mapping2);
        return;
    }

    unsigned int best_start = 0;
    unsigned int least_diff = diffs[0];
    for(int i=1; i < width-2; i++) {
        if (diffs[i] < least_diff) {
            best_start = i;
            least_diff = diffs[i];
        }
    }

    unsigned int left = best_start, right = best_start;
    unsigned int acc[4] = {row[best_start].r,row[best_start].g,row[best_start].b,row[best_start].a};

    bool left_stuck = false, right_stuck = false;
    do {
        unsigned int len = right - left + 1;
        const rgba_pixel avg = (rgba_pixel){
            .r = acc[0] / len,
            .g = acc[1] / len,
            .b = acc[2] / len,
            .a = acc[3] / len,
        };

        if (!left_stuck && left > 0) {
            rgba_pixel x0 = row[left-1];
            unsigned int diff = colordifference(avg, x0);
            if (diff < maxerror) {
                acc[0] += x0.r;
                acc[1] += x0.g;
                acc[2] += x0.b;
                acc[3] += x0.a;
                left--;
            }
            else left_stuck = true;
        }
        else left_stuck = true;

        if (!right_stuck && right < width-2) {
            rgba_pixel x0 = row[right+1];
            unsigned int diff = colordifference(avg, x0);
            if (diff < maxerror) {
                acc[0] += x0.r;
                acc[1] += x0.g;
                acc[2] += x0.b;
                acc[3] += x0.a;
                right++;
            }
            else right_stuck = true;
        }
        else right_stuck = true;

    } while(!left_stuck || !right_stuck);

    unsigned int len = right - left + 1;
    if (len > 16) {
        rgba_pixel avg = (rgba_pixel){
            .r = acc[0] / len,
            .g = acc[1] / len,
            .b = acc[2] / len,
            .a = acc[3] / len,
        };
        for(int x=left; x <= right; x++) {
            row[x] = avg;
        }
    } else {
        //left = right = best_start;
        remap_line(row+left, len, mapping1, mapping2);
    }

    if (right < width-1) {
        rle_line(row+right+1, width-right-1, diffs+right+1, maxerror, mapping1, mapping2);
    }
    if (left > 0) {
        rle_line(row, left, diffs, maxerror, mapping1, mapping2);
    }
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
    double maxerror = 0;

    int ch;
    while ((ch = getopt(argc, argv, "hvdq:Q:")) != -1) {
        switch (ch) {
            case 'd': dither = true; break;
            case 'v': verbose = true; break;
            case 'q':
            case 'Q':
                maxerror = quality_to_mse(atol(optarg));
                break;
            case '?': case 'h':
            default:
                usage(argv[0]);
                return 1;
        }
    }
    int argn = optind;

    int maxcolors = maxerror > 0 ? 255 : 0;
    if (argc==(argn+1)) {
        maxcolors=atoi(argv[argn]);
        argn++;
    }

    if (argc != argn || maxcolors < 2 || maxcolors > 255) {
        usage(argv[0]);
        return 1;
    }

    setmode(1, O_BINARY);
    setmode(0, O_BINARY);

    png24_image img;
    pngquant_error retval;

    if ((retval = rwpng_read_image24(stdin, &img))) {
        fprintf(stderr, "Error: cannot read PNG from stdin\n");
        return retval;
    }
    image_gamma = 1.0/img.gamma;

    posterize(&img, maxcolors, maxerror, dither, verbose);

    if ((retval = rwpng_write_image24(stdout, &img))) {
        fprintf(stderr, "Error: cannot write PNG to stdout\n");
        return retval;
    }

    return 0;
}


static void posterize(png24_image *img, unsigned int maxcolors, const double maxerror, bool dither, bool verbose)
{
    hist_entry histogram[256]={{0}};
    intensity_histogram(img, histogram);

    // reserve colors for black and white
    // and omit them from histogram to avoid confusing median cut
    unsigned int reservedcolors=0;
    if (BOTH(histogram[0]) >= 1.0 && maxcolors > 2) {
        maxcolors--;reservedcolors++;
        histogram[0]=(hist_entry){0,0};
    }
    if (BOTH(histogram[255]) >= 1.0 && maxcolors > 2) {
        maxcolors--;reservedcolors++;
        histogram[255]=(hist_entry){0,0};
    }

    palette pal;
    unsigned int levels = reduce(maxcolors, maxerror, histogram, &pal);

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

    double rle_treshold = (maxerror ? maxerror : last_err) * 0.8;
    const double ugly_treshold = quality_to_mse(80);
    if (rle_treshold > ugly_treshold) {
        rle_treshold = (rle_treshold + ugly_treshold)/2.0;
    }
    const double awful_treshold = quality_to_mse(50);
    if (rle_treshold > awful_treshold) {
        rle_treshold = awful_treshold;
    }

    if (verbose) {
    fprintf(stderr, "RLE at q=%d\n", mse_to_quality(rle_treshold));
    }


    rle(img, rle_treshold, &pal, dither);
}
