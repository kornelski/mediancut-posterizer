#Median cut posterizer

Reduces number of distinct color/alpha intensities in the image. Unlike typical posterization, which distributes levels evenly, this one tries to pick best levels using Median Cut.

The goal of this tool is to make RGB/RGBA images more compressible, assuming that lower number of unique byte values increses chance of finding repetition and improves efficiency of Huffman coding.

##Usage

    posterize 30 < input.png > output.png

Only stdin/stdout is supported. First (only) argument is number of levels to use, in addition to 0 and 255 which are always used.