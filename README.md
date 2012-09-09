#Median Cut PNG Posterizer

Reduces number of distinct color/alpha intensities in the image. Unlike typical posterization, which distributes levels evenly, this one tries to pick levels intelligently using varaince-based Median Cut and Voronoi iteration.

The goal of this tool is to make RGB/RGBA PNG images more compressible, assuming that lower number of unique byte values increses chance of finding repetition and improves efficiency of Huffman coding.

##Usage

    posterize [ -v ] [ -d ] [ -q <quality> ] [ levels ] < input.png > output.png

* `levels` — Number of levels to use (2-255). Lower number gives worse quality, but smaller file.
* `-q num` — Picks minimum number of levels needed to achieve given quality. `num` is quality 0-100 (100 is best, similar to JPEG). Number of levels is optional if quality is specified.
* `-d` — Enables simple ordered dithering.
* `-v` — Verbose output. Prints mean square error (MSE) caused by posterization.

Only stdin/stdout is supported.

##GUI?

Integrated in [ImageAlpha.app](http://pngmini.com).
