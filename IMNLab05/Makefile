.ONESHELL:
CC=gcc
CFLAGS= -O2 -lm -Wall
DEPS=
OBJS=main.o
OUT=main

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CFLAGS)

.c.o:
	$(CC) -c -o $@ $< $(CFLAGS)

run: $(OUT)
	./$(OUT)

all:
	make clean
	make run

plot: $(OUT)
	gnuplot gnuplot.sh

clean:
	rm -f *o
	rm -f *~
	rm -f *.dat
	rm -f *.eps
	rm -f *.png
	rm -f $(OUT)

.PHONY: clean plot all
