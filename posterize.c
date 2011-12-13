
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "png.h"
#include "rwpng.h"

static void interpolate_palette_front(int palette[], int dither);

struct box {
    double sum, variance;
    int start, end;
};

static double weighted_avg(struct box box, double histogram[])
{
    double weight=0,sum=0;
    for(int val=box.start; val < box.end; val++) {
        weight += histogram[val];
        sum += val*histogram[val];
    }
    return weight ? sum/weight : 0.0f;
}

static double variance(struct box box, double histogram[])
{
    double avg = weighted_avg(box,histogram);

    double weight=0,sum=0;
    for(int val=box.start; val < box.end; val++) {
        weight += histogram[val];
        sum += (avg-val)*(avg-val)*histogram[val];
    }
    return weight ? sum/weight : 0.0f;
}

static double palette_mse(double histogram[], int palette[])
{
    int tmp[256];
    memcpy(tmp, palette, 256*sizeof(palette[0]));
    interpolate_palette_front(tmp, 0);

    double mse=0, hist_total=0;
    for (int i=0; i < 256; i++) {
        int best = tmp[i];
        mse += (i-best)*(i-best) * histogram[i];

        hist_total += histogram[i];
    }
    return mse / hist_total;
}

static void palette_from_boxes(struct box boxes[], int numboxes, double histogram[], int palette[])
{
    memset(palette, 0, 256*sizeof(palette[0]));

    for(int box=0; box < numboxes; box++) {
        int value = round(weighted_avg(boxes[box],histogram));
        palette[value] = value;
    }
}

/*
 1-dimensional median cut, using variance for "largest" box
*/
static void reduce(const int maxcolors, double histogram[], int palette[])
{
    int numboxes=1;
    struct box boxes[256];

    boxes[0].start=1; // skip first and last entry, as they're always included
    boxes[0].end=255;
    boxes[0].sum=0;
    for(int i=boxes[0].start; i < boxes[0].end; i++) boxes[0].sum += histogram[i];
    boxes[0].variance = variance(boxes[0], histogram);

    while(numboxes < maxcolors) {
        int boxtosplit=-1;
        int largest=0;
        for(int box=0; box < numboxes; box++) {
            if (boxes[box].variance > largest && (boxes[box].end-boxes[box].start)>=2) {
                largest = boxes[box].variance;
                boxtosplit=box;
            }
        }
        if (boxtosplit < 0) {
            break;
        }
        double sum=0;
        int val=boxes[boxtosplit].start;
        for(; val < boxes[boxtosplit].end-1; val++) {
            sum += histogram[val];
            if (sum >= boxes[boxtosplit].sum/2.0) {
                break;
            }
        }

        boxes[numboxes].start = boxes[boxtosplit].start;
        boxes[numboxes].end = val+1;
        boxes[numboxes].sum = sum;
        boxes[numboxes].variance = variance(boxes[numboxes], histogram);
        boxes[boxtosplit].start = val+1;
        boxes[boxtosplit].sum -= boxes[numboxes].sum;
        boxes[boxtosplit].variance = variance(boxes[boxtosplit], histogram);
        numboxes++;
    }

    palette_from_boxes(boxes, numboxes, histogram, palette);
}

static void remap(read_info img, const int *palette1, const int *palette2)
{
    for(int i=0; i < img.height; i++) {
        for(int j=0; j < img.width; j++) {
            int x = j*4;
            const int *palette = (i^j)&1 ? palette1 : palette2;

            int a = palette[img.row_pointers[i][x+3]];
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

static void intensity_histogram(read_info img, double histogram[])
{
    for(int i=0; i < img.height; i++) {
        for(int x=0; x < img.width*4; x+=4) {
            double a = img.row_pointers[i][x+3]/255.0;

            // opaque colors get more weight
            histogram[img.row_pointers[i][x]] += a;
            histogram[img.row_pointers[i][x+1]] += a;
            histogram[img.row_pointers[i][x+2]] += a;
            histogram[img.row_pointers[i][x+3]] += 1.0 - a;
        }
    }
}

static void interpolate_palette_front(int palette[], int dither)
{
    int nextval=0, lastval=0;
    palette[0]=0;
    palette[255]=255; // 0 and 255 are always included

    for(int val=0; val < 256; val++) {
        if (palette[val]==val) {
            lastval = val;
            for(int j=val+1; j < 256; j++) {
                if (palette[j]==j) {nextval=j; break;}
            }
        }
        if (!dither) {
            palette[val] = (val - lastval) < (nextval - val) ? lastval : nextval;
        } else {
            palette[val] = (val - lastval)/2 < (nextval - val) ? lastval : nextval;
        }
    }
}

static void interpolate_palette_back(int palette2[], int dither)
{
    int nextval=255, lastval=255;
    palette2[0]=0;
    palette2[255]=255; // 0 and 255 are always included

    for(int val=255; val >=0; val--) {
        if (palette2[val]==val) {
            lastval = val;
            for(int j=val-1; j >= 0; j--) {
                if (palette2[j]==j) {nextval=j; break;}
            }
        }
        if (!dither) {
            palette2[val] = (val - lastval) >= (nextval - val) ? lastval : nextval;
        } else {
            palette2[val] = (val - lastval)/2 >= (nextval - val) ? lastval : nextval;
        }
    }
}

static void dither_palette(int palette[], int palette2[], int dither)
{
    memcpy(palette2, palette, sizeof(int)*256);

    // front to back. When dithering, it's biased towards nextval
    interpolate_palette_front(palette, dither);

    // back to front, so dithering bias is the other way.
    interpolate_palette_back(palette2, dither);
}

static void usage(const char *exepath)
{
    const char *name = strrchr(exepath, '/');
    if (name) name++; else name = exepath;
    fprintf(stderr, "Median Cut PNG Posterizer 1.2 (2011).\n" \
            "Usage: %s [-d] levels\n\n" \
            "Specify number of levels 2-255 as an argument. -d enables dithering\n" \
            "Image is always read from stdin and written to stdout.\n"
            "%s -d 16 < in.png > out.png\n", name, name);
}

static void voronoi(double histogram[], int palette[])
{
    interpolate_palette_front(palette, 0);

    double counts[256] = {0};
    double sums[256] = {0};

    // remap palette
    for (int i=0; i < 256; i++) {
        int best = palette[i];
        counts[best] += histogram[i];
        sums[best] += histogram[i] * (double)i;
    }

    memset(palette, 0, 256*sizeof(palette[0]));

    // rebuild palette from remapped averages
    for(int i=0; i < 256; i++) {
        if (counts[i]) {
            int value = round(sums[i]/counts[i]);
            palette[value] = value;
        }
    }
}

int main(int argc, char *argv[])
{
    int argn=1;
    int dither=0;
    if (argc==3 && 0==strcmp("-d", argv[1])) {
        dither=1;
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

    int palette[256], palette2[256];
    reduce(maxcolors, histogram, palette);

    double last_mse = INFINITY;
    for(int j=0; j < 100; j++) {
        voronoi(histogram, palette);
        double new_mse = palette_mse(histogram, palette);
        if (new_mse == last_mse) break;
        last_mse = new_mse;
    }

    dither_palette(palette, palette2, dither);

    remap(img, palette, palette2);

    if ((retval = rwpng_write_image_init(stdout, &img)) ||
        (retval = rwpng_write_image_whole(&img))) {
        fprintf(stderr, "Error: cannot write PNG to stdout\n");
        return retval;
    }

    return 0;
}
