#Median Cut PNG Posterizer

Reduces number of distinct color/alpha intensities in the image. Unlike typical posterization, which distributes levels evenly, this one tries to pick levels intelligently using varaince-based Median Cut and Voronoi iteration.

The goal of this tool is to make RGB/RGBA PNG images more compressible, assuming that lower number of unique byte values increses chance of finding repetition and improves efficiency of Huffman coding.

##Usage

    posterize [ -v ] [ -d ] [ -Q <quality> ] [ levels ] -- [ input.png ] [ output.png ]

* `levels` — Number of levels to use (2-255). Lower number gives worse quality, but smaller file.
* `-Q num` — Picks minimum number of levels needed to achieve given quality. `num` is quality 0-100 (100 is best, similar to JPEG). Number of levels is optional if quality is specified.
* `-d` — Enables simple ordered dithering.
* `-v` — Verbose output. Prints mean square error (MSE) caused by posterization.

If input/output files are not specified then stdin/stdout is used respectively.

Posterized images can be further compressed using [PNGOUT](http://www.jonof.id.au/kenutils) or similar. Try [ImageOptim](http://imageoptim.com).

##GUI?

Integrated in [ImageAlpha.app](http://pngmini.com).

##Licenses

###Posterizer

© 2011-2012 Kornel Lesiński.

This program is free software: you can redistribute it and/or modify
it under the terms of the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html)
as published by the Free Software Foundation, either version 3
of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


###`rwpng.c`

© 1998-2000 Greg Roelofs.  All rights reserved.

This software is provided "as is," without warranty of any kind,
express or implied.  In no event shall the author or contributors
be held liable for any damages arising in any way from the use of
this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute
it freely, subject to the following restrictions:

1. Redistributions of source code must retain the above copyright
 notice, disclaimer, and this list of conditions.
2. Redistributions in binary form must reproduce the above copyright
 notice, disclaimer, and this list of conditions in the documenta-
 tion and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this
 software must display the following acknowledgment:

   > This product includes software developed by Greg Roelofs
   > and contributors for the book, "PNG: The Definitive Guide,"
   > published by O'Reilly and Associates.
