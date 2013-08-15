CFLAGSOPT ?= -DNDEBUG -O3 -fstrict-aliasing -ffast-math -funroll-loops -fomit-frame-pointer -ffinite-math-only

CFLAGS ?= -Wall -Wno-unknown-pragmas -I. -I/usr/local/include/ -I/usr/include/ -I/usr/X11/include/ $(CFLAGSOPT)
CFLAGS += -std=c99 $(CFLAGSADD)

LDFLAGS ?= -L/usr/local/lib/ -L/usr/lib/ -L/usr/X11/lib/
LDFLAGS += -lpng -lz -lm $(LDFLAGSADD)

OBJS=rwpng.o
COCOA_OBJS = rwpng_cocoa.o

ifdef USE_COCOA
CFLAGS += -DUSE_COCOA=1
OBJS += $(COCOA_OBJS)
FRAMEWORKS += -framework Cocoa
endif

all: posterize blurize

blurize: blurize.o $(OBJS)
	$(CC) blurize.o $(OBJS) $(LDFLAGS) $(FRAMEWORKS) -o $@

posterize: posterize.o $(OBJS)
	$(CC) posterize.o $(OBJS) $(LDFLAGS) $(FRAMEWORKS) -o $@

rwpng_cocoa.o: rwpng_cocoa.m
	clang -c $(CFLAGS) -o $@ $<

clean:
	-rm -f $(OBJS)
