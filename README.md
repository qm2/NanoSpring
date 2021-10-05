# NanoSpring

[![Build Status](https://travis-ci.com/qm2/NanoSpring.svg?branch=master)](https://travis-ci.com/qm2/NanoSpring)

NanoSpring - Tool for compression of nanopore genomic reads in FASTQ format (gzipped input also supported). Compresses only the read sequences (i.e., ignores quality values and read identifiers). Achieves over 3 times better compression ratios than Gzip (for recent basecaller versions) and fast decompression. The algorithm requires C++11 and g++ compiler and works on Linux (tested on Ubuntu 16.04, 18.04 and 20.04). 

#### BioRxiv: https://www.biorxiv.org/content/10.1101/2021.06.09.447198v1

## Install from source code

### Download repository

```
git clone --recursive https://github.com/qm2/NanoSpring.git
```

### Install
The instructions below will create the NanoSpring executable in the build directory inside NanoSpring. If you plan to build and run NanoSpring on an older architectures, then you might need to remove the flag ```-msse4.1``` from the ```target_compile_options``` in CMakeLists.txt (or use flags based on the target architecture).

On Linux with cmake installed and version at least 3.12 (check using `cmake --version`):
```
cd NanoSpring
mkdir build
cd build
cmake ..
make -j
```
On Linux with cmake not installed or with version older than 3.12:
```bash
cd NanoSpring
mkdir build
cd build
wget https://cmake.org/files/v3.12/cmake-3.12.4.tar.gz
tar -xzf cmake-3.12.4.tar.gz
cd cmake-3.12.4
./configure
make -j
cd ..
./cmake-3.12.4/bin/cmake ..
make -j
```

Test installation by running the test script (from the `build` directory):
```
./../util/test_script.sh
```

## Usage
Run the NanoSpring compressor executable with the options below:
```
Allowed options:
  -h [ --help ]                  produce help message
  -c [ --compress ]              compress
  -d [ --decompress ]            decompress
  -i [ --input-file ] arg        input file name
  -o [ --output-file ] arg       output file name
  -t [ --num-threads ] arg (=20) number of threads (default 20)
  -k [ --kmer ] arg (=23)        kmer size for the minhash (default 23)
  -n [ --num-hash ] arg (=60)    number of hash functions for minhash (default
                                 60)
  --overlap-sketch-thr arg (=6)  the overlap sketch threshold for minhash
                                 (default 6)
  --minimap-k arg (=20)          kmer size for the minimap2 (default 20)
  --minimap-w arg (=50)          window size for the minimap2 (default 50)
  --max-chain-iter arg (=400)    the max number of partial chains during
                                 chaining for minimap2 (default 400)
  --edge-thr arg (=4000000)      the max number of edges allowed in a consensus
                                 graph (default 4000000)
  -w [ --working-dir ] arg (=.)  directory to create temporary files (default
                                 current directory)
```
Note that the compressed files are tar archives consisting of the different compressed streams, although we recommend using the ```.NanoSpring``` extension as in the examples shown below.


## Example Usage of NanoSpring
This section contains several examples for compression and decompression with various modes and options. The compressed file uses the ```.NanoSpring``` extension as a convention. The decompressed file uses the ```.reads``` extension as a convention.

For compressing file.fastq using default 20 threads.
```
./NanoSpring -c -i file.fastq -o file.NanoSpring
```
Using 10 threads.
```
./NanoSpring -c -i file.fastq -o file.NanoSpring -t 10
```
For compressing file.fastq.gz (gzipped fastq files) using default 20 threads.
```
./NanoSpring -c -i file.fastq.gz -o file.NanoSpring 
```
Decompressing the file.NanoSpring with default 20 threads to file.reads.
```
./NanoSpring -d -i file.NanoSpring -o file.reads
```
Compare decompressed file to reads in original file to verify compression was lossless.
```
cmp file.reads <(cat file.fastq | sed -n '2~4p')
```
Decompressing the file.NanoSpring with 10 threads to file.reads.
```
./NanoSpring -d -i file.NanoSpring -o file.reads -t 10
```
Compressing with kmer size for minhash to be 20, number of hash functions for minhash to be 64, and the overlap sketch threhold for minhash to be 8.
```
./NanoSpring -c -i file.fastq -k 20 -n 64 --overlap-sketch-thr 8 -o file.NanoSpring 
```
Compressing with kmer size of minimap2 to be 25 and the window size for the minimap2 to be 100.
```
./NanoSpring -c -i file.fastq --minimap-k 25 --minimap-w 100 -o file.NanoSpring
```
Compressing with the max number of partial chains during chaining for minimap2 to be 800.
```
./NanoSpring -c -i file.fastq --max-chain-iter 800 -o file.NanoSpring  
```
Compressing with the max number of edges allowed in a consensus graph to be 8000000.
```
./NanoSpring -c -i file.fastq --edge-thr 8000000 -o file.NanoSpring  
```
