# NanoSpring

## Install from source code

### Download repository

```
git --recursive clone https://github.com/fanzhuyifan/NanoSpring.git
```

### Install
The instructions below will create the executables (both compressor and decompressor) in the build directory inside NanoSpring.

On Linux with cmake installed and version at least 3.10 (check using cmake --version):
```
cd NanoSpring
mkdir build
cd build
cmake ..
make
```
## Usage
### Compression - compresses FASTQ reads. 

Run the NanoSpring compressor executable with the options below:
```
Allowed options:
  -h [ --help ]                         produce help message
  -i [ --input-file ] arg               input file name
  -o [ --output-file ] arg (=compressedFile)
                                        output file name
  -t [ --num-threads ] arg (=16)        number of threads (default 16)
  -k [ --kmer ] arg (=23)               kmer size for the minhash (default 23)
  -n [ --num-hash ] arg (=64)           number of hash functions for minhash
                                        (default 64)
  --overlap-sketch-thr arg (=8)         the overlap sketch threhold for minhash
                                        (default 8)
  --minimap-k arg (=20)                 kmer size for the minimap2 (default 20)
  --minimap-w arg (=50)                 window size for the minimap2 (default
                                        50)
  --max-chain-iter arg (=400)           the max number of partial chains during
                                        chaining for minimap2 (default 400)
  --edge-thr arg (=2000000)             the max number of edges allowed in a
                                        consensus graph
  -g [ --gzipped-fastq ]                enable if compression input is gzipped
                                        fastq
  -w [ --working-dir ] arg (=.)         directory to create temporary files
                                        (default current directory)
```
Note that the compressed files are tar archives consisting of the different compressed streams, although we recommend using the .NanoSpring extension as in the examples shown below.

### Decompression -  decompresses reads. 

Run the NanoSpring decompressor executable with the options below:
```
Allowed options:
  -h [ --help ]                         produce help message
  -i [ --input-file ] arg               input file name
  -o [ --output-file ] arg (=compressedFile)
                                        output file name
  -t [ --num-threads ] arg (=16)        number of threads (default 16)
```

## Example Usage of Nanopore Compression
This section contains several examples for compression and decompression with various modes and options. The compressed file uses the .NanoSpring extension as a convention.

For compressing file.fastq losslessly using default 16 threads (Lossless).
```
./testCompressor -i file.fastq -o file.NanoSpring
```
Using 20 threads (Lossless).
```
./testCompressor -i file.fastq -o file.NanoSpring -t 20
```
For compressing file.fastq.gz (gzipped fastq files) losslessly using default 16 threads (Lossless).
```
./testCompressor -i file.fastq -o file.NanoSpring -g 
```
Compressing with kmer size for minhash 23, number of hash functions for minhash 54, and the overlap sketch threhold for minhash 6.
```
./testCompressor -i file.fastq -k 23 -n 54 --overlap-sketch-thr 6 -o file.NanoSpring 
```
Compressing with kmer size of minimap2 25 and the window size for the minimap2 100.
```
./testCompressor -i file.fastq --minimap-k 25 --minimap-w 100 -o file.NanoSpring
```
Compressing with the max number of partial chains during chaining for minimap2 to be 800.
```
./testCompressor -i file.fastq --max-chain-iter 800 -o file.NanoSpring  
```
Compressing with the max number of edges allowed in a consensus graph to be 4000000.
```
./testCompressor -i file.fastq --edge-thr 4000000 -o file.NanoSpring  
```
Decompressing with default 16 threads to file.reads.
```
./testDecompressor -i file.NanoSpring -o file.reads
```
