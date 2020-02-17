
OS_NAME := $(shell uname -s | tr A-Z a-z)

ifeq ($(OS_NAME),linux)
	CC := g++
else
	CC := g++-9
endif

all: enano

enano: 
		$(CC) -fopenmp enano_fastq.cpp Compressor.cpp -o enano

clean:
		rm -f enano *.o
