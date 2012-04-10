
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "png.h"
#include "rwpng.h"

static void interpolate_palette_front(unsigned int palette[], const bool dither);

// Converts gamma 2.0 (approx of 2.2) to linear unit value. Linear color is required for preserving brightness (esp. when dithering).
inline static double gamma_to_linear(const double value)
{
    return sqrt(value/255.0);
}

// Reverses gamma_to_linear. *256 is not off-by-one error.
inline static unsigned int linear_to_gamma(const double value)
{
    const double g = value*value*256.0;
    return g < 255.0 ? g : 255;
}


// median cut "box" in this implementation is actually a line,
// since it only needs to track lowest/highest intensity
struct box {
    double sum, variance;
    unsigned int start, end;
};

// average values in a "box" proportionally to frequency of their occurence
static double weighted_avg_linear(const unsigned int start, const unsigned int end, const double histogram[])
{
    double weight=0,sum=0;
    for(unsigned int val=start; val < end; val++) {
        weight += histogram[val];
        sum += gamma_to_linear(val)*histogram[val];
    }
    return weight ? sum/weight : 0;
}

// variance (AKA second moment) of the box. Measures how much "spread" the values are
static double variance_in_range(const unsigned int start, const unsigned int end, const double histogram[])
{
    const double avg = weighted_avg_linear(start, end, histogram);

    double sum=0;
    for(unsigned int val=start; val < end; val++) {
        const double delta = avg-gamma_to_linear(val);
        sum += delta*delta*histogram[val];
    }
    return sum;
}

static double variance(const struct box box, const double histogram[])
{
    return variance_in_range(box.start, box.end, histogram);
}

// Square error. Estimates how well palette "fits" the histogram.
static double palette_error(const double histogram[], const unsigned int palette_orig[])
{
    unsigned int palette[256];
    memcpy(palette, palette_orig, 256*sizeof(palette[0]));

    // the input palette has gaps
    interpolate_palette_front(palette, false);

    double se=0;
    for (int i=0; i < 256; i++) {
        double delta = gamma_to_linear(i)-gamma_to_linear(palette[i]);
        se += delta*delta*histogram[i];
    }
    return se;
}

// converts boxes to palette.
// palette here is a sparse array where elem[x]=x is taken, elem[x]=0 is free (except x=0)
static void palette_from_boxes(const struct box boxes[], const int numboxes, const double histogram[], unsigned int palette[])
{
    memset(palette, 0, 256*sizeof(palette[0]));

    for(int box=0; box < numboxes; box++) {
        int value = linear_to_gamma(weighted_avg_linear(boxes[box].start, boxes[box].end, histogram));
        palette[value] = value;
    }
}

/*
 1-dimensional median cut, using variance for "largest" box
*/
static void reduce(const unsigned int maxcolors, const double histogram[], unsigned int palette[])
{
    unsigned int numboxes=1;
    struct box boxes[256];

    // build the first "box" that encompasses all values
    boxes[0].start=1; // skip first and last entry, as they're always included
    boxes[0].end=255;
    boxes[0].sum=0;
    for(unsigned int i=boxes[0].start; i < boxes[0].end; i++) boxes[0].sum += histogram[i];
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
        for(unsigned int i=boxes[boxtosplit].start; i < bestsplit; i++) sum += histogram[i];

        // create new boxes from halves
        boxes[numboxes].start = boxes[boxtosplit].start;
        boxes[numboxes].end = bestsplit;
        boxes[numboxes].sum = sum;
        boxes[numboxes].variance = variance(boxes[numboxes], histogram);
        boxes[boxtosplit].start = bestsplit;
        boxes[boxtosplit].sum -= boxes[numboxes].sum;
        boxes[boxtosplit].variance = variance(boxes[boxtosplit], histogram);
        numboxes++;
    }

    palette_from_boxes(boxes, numboxes, histogram, palette);
}

// palette1/2 is for even/odd pixels, allowing very simple "ordered" dithering
static void remap(read_info img, const unsigned int palette1[], const unsigned int palette2[])
{
    for(unsigned int i=0; i < img.height; i++) {
        for(unsigned int j=0; j < img.width; j++) {
            unsigned int x = j*4;
            const unsigned int *palette = (i^j)&1 ? palette1 : palette2;

            const unsigned int a = palette[img.row_pointers[i][x+3]];
            if (a) {
                img.row_pointers[i][x] = palette[img.row_pointers[i][x]];
                img.row_pointers[i][x+1] = palette[img.row_pointers[i][x+1]];
                img.row_pointers[i][x+2] = palette[img.row_pointers[i][x+2]];
                img.row_pointers[i][x+3] = a;
            } else {
                // clear "dirty alpha"
                img.row_pointers[i][x] = 0;
                img.row_pointers[i][x+1] = 0;
                img.row_pointers[i][x+2] = 0;
                img.row_pointers[i][x+3] = 0;
            }
        }
    }
}

// usually RGBA images are stored/rendered in "premultiplied" format which is R*A, G*A, B*A
// this causes loss of precision, so it may be a good idea to posterize to this value anyway
inline static unsigned int premultiplied_alpha_rounding(const unsigned int value, const unsigned int alpha)
{
    return value * alpha / alpha;
}

// it doesn't count unique colors, only intensity values of all channels
static void intensity_histogram(const read_info img, double histogram[])
{
    for(unsigned int i=0; i < img.height; i++) {
        const unsigned char *const row = img.row_pointers[i];
        for(unsigned int x=0; x < img.width*4; x+=4) {
            const unsigned int alpha = row[x+3];
            if (alpha) {
                // opaque colors get more weight
                const double weight = alpha/255.0;

                histogram[premultiplied_alpha_rounding(row[x], alpha)] += weight;
                histogram[premultiplied_alpha_rounding(row[x+1], alpha)] += weight;
                histogram[premultiplied_alpha_rounding(row[x+2], alpha)] += weight;
                histogram[row[x+3]] += 1.0;
            }
            else histogram[0] += 4.0;
        }
    }
}

// interpolates front-to-back. If dither is true, it will bias towards one side
static void interpolate_palette_front(unsigned int palette[], const bool dither)
{
    unsigned int nextval=0, lastval=0;
    palette[0]=0;
    palette[255]=255; // 0 and 255 are always included

    for(unsigned int val=0; val < 256; val++) {
        if (palette[val]==val) {
            lastval = val;
            for(unsigned int j=val+1; j < 256; j++) {
                if (palette[j]==j) {nextval=j; break;}
            }
        }
        const double lastvaldiff = (gamma_to_linear(val) - gamma_to_linear(lastval));
        const double nextvaldiff = (gamma_to_linear(nextval) - gamma_to_linear(val));
        if (!dither) {
            palette[val] = lastvaldiff < nextvaldiff ? lastval : nextval;
        } else {
            palette[val] = lastvaldiff/2 < nextvaldiff ? lastval : nextval;
        }
    }
}

// interpolates back-to-front. Always biased for dither.
static void interpolate_palette_back(unsigned int palette2[])
{
    unsigned int nextval=255, lastval=255;
    palette2[0]=0;
    palette2[255]=255; // 0 and 255 are always included

    for(int val=255; val >= 0; val--) {
        if (palette2[val]==val) {
            lastval = val;
            for(int j=val-1; j >= 0; j--) {
                if (palette2[j]==j) {nextval=j; break;}
            }
        }
        const double lastvaldiff = (gamma_to_linear(val) - gamma_to_linear(lastval));
        const double nextvaldiff = (gamma_to_linear(nextval) - gamma_to_linear(val));
        palette2[val] = lastvaldiff/2 >= nextvaldiff ? lastval : nextval;
    }
}

static void dither_palette(unsigned int palette[], unsigned int palette2[])
{
    memcpy(palette2, palette, sizeof(unsigned int)*256);

    // front to back. When dithering, it's biased towards nextval
    interpolate_palette_front(palette, true);

    // back to front, so dithering bias is the other way.
    interpolate_palette_back(palette2);
}

static void usage(const char *exepath)
{
    const char *name = strrchr(exepath, '/');
    if (name) name++; else name = exepath;
    fprintf(stderr, "Median Cut PNG Posterizer 1.3 (2012).\n" \
            "Usage: %s [-d] levels\n\n" \
            "Specify number of levels 2-255 as an argument. -d enables dithering\n" \
            "Image is always read from stdin and written to stdout.\n"
            "%s -d 16 < in.png > out.png\n", name, name);
}

// performs voronoi iteration (mapping histogram to palette and creating new palette from remapped values)
// this shifts palette towards local optimum
static void voronoi(const double histogram[], unsigned int palette[])
{
    interpolate_palette_front(palette, false);

    double counts[256] = {0};
    double sums[256] = {0};

    // remap palette
    for (unsigned int i=0; i < 256; i++) {
        int best = palette[i];
        counts[best] += histogram[i];
        sums[best] += histogram[i] * (double)i;
    }

    memset(palette, 0, 256*sizeof(palette[0]));

    // rebuild palette from remapped averages
    for(unsigned int i=0; i < 256; i++) {
        if (counts[i]) {
            int value = round(sums[i]/counts[i]);
            palette[value] = value;
        }
    }
}

int main(int argc, char *argv[])
{
    int argn=1;
    bool dither = false;
    if (argc==3 && 0==strcmp("-d", argv[1])) {
        dither=true;
        argn++;
    }
    int maxcolors=0;
    if (argc==(argn+1)) {
        maxcolors=atoi(argv[argn]);
        argn++;
    }

    if (argc != argn || maxcolors < 2 || maxcolors > 255) {
        usage(argv[0]);
        return 1;
    }

    read_info img;
    pngquant_error retval;
    if ((retval = rwpng_read_image(stdin, &img))) {
        fprintf(stderr, "Error: cannot read PNG from stdin\n");
        return retval;
    }

    double histogram[256]={0};
    intensity_histogram(img, histogram);

    // reserve colors for black and white
    // and omit them from histogram to avoid confusing median cut
    if (histogram[0] && maxcolors>2) {maxcolors--; histogram[0]=0;}
    if (histogram[255] && maxcolors>2) {maxcolors--; histogram[255]=0;}

    unsigned int palette[256], palette2[256];
    reduce(maxcolors, histogram, palette);

    double last_err = INFINITY;
    for(unsigned int j=0; j < 100; j++) {
        voronoi(histogram, palette);
        double new_err = palette_error(histogram, palette);
        if (new_err == last_err) break;
        last_err = new_err;
    }

    if (dither) {
        dither_palette(palette, palette2);
        remap(img, palette, palette2);
    } else {
        interpolate_palette_front(palette, false);
        remap(img, palette, palette);
    }

    if ((retval = rwpng_write_image_init(stdout, &img)) ||
        (retval = rwpng_write_image_whole(&img))) {
        fprintf(stderr, "Error: cannot write PNG to stdout\n");
        return retval;
    }

    return 0;
}
