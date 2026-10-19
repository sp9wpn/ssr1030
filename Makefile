CC      ?= gcc
CFLAGS  ?= -O2 -g -Wall -Wextra -std=gnu11
HACKRF_CFLAGS := $(shell pkg-config --cflags libhackrf 2>/dev/null | sed 's|/libhackrf||')
HACKRF_LIBS   := $(shell pkg-config --libs libhackrf 2>/dev/null || echo -lhackrf)
LDLIBS  = $(HACKRF_LIBS) -lpthread -lm

SRC = src/main.c src/demod.c src/modes.c src/analysis.c src/display.c src/output.c src/adsb.c
OBJ = $(SRC:.c=.o)

ssr1030: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDLIBS)

src/%.o: src/%.c src/ssr1030.h src/analysis.h
	$(CC) $(CFLAGS) $(HACKRF_CFLAGS) -c -o $@ $<

clean:
	rm -f ssr1030 src/*.o

.PHONY: clean
