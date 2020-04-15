
OS_NAME := $(shell uname -s | tr A-Z a-z)

OP_FLAGS := -std=c++11 -O3 -fstrict-aliasing -ffast-math -fomit-frame-pointer

ifeq ($(OS_NAME),linux)
	CC := g++
else
	CC := g++-9
endif

all: enano

enano: 
		$(CC) $(OP_FLAGS) -fopenmp enano_fastq.cpp Compressor.cpp -o enano

clean:
		rm -f enano *.o
