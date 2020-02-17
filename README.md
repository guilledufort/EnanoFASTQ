# ENANO FASTQ
## An encoder for nanopore FASTQ files
#### Publication: on review
## Description
ENANO is a novel FASTQ compression algorithm especially designed for nanopore sequencing FASTQ files. We tested ENANO and current state-of-the-art compressors on several publicly available nanopore datasets. The results show that our algorithm consistently achieves the best compression performance on every nanopore dataset, while being computationally efficient in terms of speed and memory requirements when compared to existing alternatives.

## Requirements
0. g++ ( >= 4.8)
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

To compile enano run:
```bash
cd EnanoFASTQ/enano
make
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

| Dataset | # of files | size (GB) | Description | Link |
|------|------|------|------|------|
 *sor\** | 4 | 124.071 | Sorghum bicolor Tx430 | https://www.nature.com/articles/s41467-018-07271-1#data-availability |
 *bra\** | 18 | 43.014 | Doubled haploid canola (Brassica napus L.) | https://www.nature.com/articles/s41598-019-45131-0#data-availability |
 *lun* | 13 | 15.239 | Human lung bacterial  metagenomic | https://www.nature.com/articles/s41587-019-0156-5#data-availability |
 *joi* | 9 | 4.672 | Infected orthopaedic devices metagenomic | https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5094-y |
 *vir\** | 10 | 4.375 | Direct RNA sequencing (HSV-1) | https://www.nature.com/articles/s41467-019-08734-9#data-availability |
 *hs1* | 1 | 249.791 | Human GM12878 Utah/Ceph cell line | https://github.com/nanopore-wgs-consortium/NA12878 |
 *hs2^* | 50 | 193.920 | Human GM12878 Utah/Ceph cell line | https://www.nature.com/articles/s41467-019-09637-5#data-availability|
 *npd\** | 336 | 113.440 | Multiple organisms | https://github.com/guidufort/DualFqz |

\*Datasets that require the SRA toolkit to be downloaded. 

^We only used the first 50 files of the dataset.

### Downloading the datasets

To download a dataset you have to run the *download_script.sh* of the specific dataset.
For example, to download *sor* run:
```bash
cd EnanoFASTQ
dataset/sor/download_script.sh
```

The scripts use the command *wget* to perform the download. 
To install *wget* on macOS run:
 ```bash
brew install wget
```
To install *wget* on Ubuntu or CentOS run:
 ```bash
sudo apt-get install wget
```

Some datasets require the SRA toolkit (2.9.6-1 release) to be downloaded. To install the SRA toolkit you can follow the instructions here https://ncbi.github.io/sra-tools/install_config.html, and place the toolkit's root-folder under the EnanoFASTQ directory, or you can run one of the scripts we built. There is a different script for each OS, so you have to choose the one corresponding to your OS.
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
To decompress with 8 threads the example compressed file:
```bash
cd EnanoFASTQ
enano/enano -d -t 8 example/SAMPLE.enano example/SAMPLE_dec.fastq
```

#### Check if decoding is successful
The output has to be empty.
```bash
cmp example/SAMPLE.fastq example/SAMPLE_dec.fastq
```
### Credits
The methods used for encoding the reads names, model frequency counters, and to do the reads parsing, are the ones proposed by James Bonefield in FQZComp, with some modifications. The range coder is derived from Eugene Shelwien.
