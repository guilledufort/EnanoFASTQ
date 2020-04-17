OS_NAME := $(shell uname -s | tr A-Z a-z)

CXXFLAGS = -fopenmp -std=c++11 -O3 -march=native -fomit-frame-pointer -fstrict-aliasing -ffast-math

ifeq ($(OS_NAME),linux)
	CXX = g++
else
	CXX = g++-9
endif

all: enano

enano: 
		$(CXX) $(CXXFLAGS) enano_fastq.cpp Compressor.cpp -o enano

clean:
		rm -f enano *.o
