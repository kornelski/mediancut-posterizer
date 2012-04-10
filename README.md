#Median Cut PNG Posterizer

Reduces number of distinct color/alpha intensities in the image. Unlike typical posterization, which distributes levels evenly, this one tries to pick levels intelligently using varaince-based Median Cut and Voronoi iteration.

The goal of this tool is to make RGB/RGBA PNG images more compressible, assuming that lower number of unique byte values increses chance of finding repetition and improves efficiency of Huffman coding.

##Usage

    posterize [ -d ] levels < input.png > output.png

* `-d` — whether to use simple ordered dithering.
* `levels` — number of levels to use (2-255). Lower number gives worse quality, but smaller file.

Only stdin/stdout is supported.

##GUI?

Integrated in [ImageAlpha.app](http://pngmini.com).
