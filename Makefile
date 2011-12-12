CFLAGS ?= -O3 -Wall

CFLAGS += -std=c99
LDFLAGS += -lpng -lz

OBJS=posterize.o rwpng.o

all: posterize

posterize: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $@

clean:
	rm -f $(OBJS)
