# ENANO FASTQ
## An encoder for nanopore FASTQ files
#### Publication: on review
## Description
ENANO is a novel FASTQ compression algorithm especially designed for nanopore sequencing FASTQ files. We tested ENANO and current state-of-the-art compressors on several publicly available nanopore datasets. The results show that our algorithm consistently achieves the best compression performance on every nanopore dataset, while being computationally efficient in terms of speed and memory requirements when compared to existing alternatives.

## Requirements
0. g++ 
1. OpenMP library

## Download

```bash
git clone https://github.com/guilledufort/EnanoFASTQ.git
```

## Install 

The following instructions will create the enano executable in the directory *enano*.
To compile enano you need to have the g++ compiler and the OpenMP library for multithreading. 

On Linux (Ubuntu or CentOS) g++ usually comes installed by default, but if not run the following:
```bash
sudo apt update
sudo apt-get install g++
```

On macOS, install GCC compiler since Clang has issues with OpenMP library:
- Install HomeBrew (https://brew.sh/)
- Install GCC (this step will be faster if Xcode command line tools are already installed using ```xcode-select --install```):
```bash
brew update
brew install gcc@9
```

- Set environment variables:
```bash
export CC=gcc-9
export CXX=g++-9
```

To compile enano run:
```bash
cd EnanoFASTQ/enano
g++ -fopenmp enano_fastq.cpp Compressor.cpp -o enano
```

# USAGE

```bash 
To compress:
  enano [options] [input_file [output_file]]

    -c             To use MAX COMPRESION MODE. Default is FAST MODE.

    -s <length>    Base sequence context length. Default is 7.

    -l <lenght>    Length of the DNA sequence context. Default is 6.

    -t <num>       Maximum number of threads allowed to use by the compressor. Default is 8.

To decompress:
   enano -d [options] foo.enano foo.fastq
    -t <num>       Maximum number of threads allowed to use by the decompressor. Default is 8.
```


## Datasets information

To test our compressor we ran experiments on the following datasets. The full information of the datasets is on our publication.

| Dataset | # of files | size (GB) | Description |
|------|------|------|------|
 *sor\** | 4 | 124.071 | Sorghum bicolor Tx430 |
 *bra\** | 4 | 43.014 | Doubled haploid canola (Brassica napus L.) |
 *lun* | 13 | 15.239 | Human lung bacterial  metagenomic |
 *joi* | 9 | 4.672 | Infected orthopaedic devices metagenomic |
 *vir\** | 10 | 4.375 | Direct RNA sequencing (HSV-1) |
 *hs1* | 1 | 249.791 | Human GM12878 Utah/Ceph cell line |
 *hs2* | 50 | 193.920 | Human GM12878 Utah/Ceph cell line |
 *npd* | 336 | 113.440 | Multiple organisms |

\*Datasets that require the SRA toolkit to be downloaded. 

### Downloading the datasets

To download a dataset you have to run the *download_script.sh* of the specific dataset.
For example, to download *sor* run:
```bash
cd EnanoFASTQ
dataset/sor/download_script.sh
```

Some datasets require for the SRA toolkit to be downloaded, and for the toolkit folder to be copied to the EnanoFASTQ directory. To install the SRA toolkit you can follow the instructions here https://ncbi.github.io/sra-tools/install_config.html, or you can run one of the scripts we built. There is a different script for each OS, so you have to choose the one corresponding to your OS.
For example, to install the SRA toolkit on macOS you can run:
 ```bash
cd EnanoFASTQ
./install_SRA_mac.sh
```

### Examples

#### Compress using ENANO
To run the compressor with 4 threads on the example file:
```bash
cd EnanoFASTQ
enano/enano -s 8 -l 5 -t 4 example/SAMPLE.fastq example/SAMPLE.enano
```
#### Decompress using ENANO
To decompress with 8 threads on the example file:
```bash
cd EnanoFASTQ
enano/enano -d -t 8 example/SAMPLE.enano example/SAMPLE_dec.fastq
```

#### Check if decoding is successful
The output has to be empty.
```bash
cmp example/SAMPLE.fastq example/SAMPLE_dec.fastq
```
