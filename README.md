## Install from source code

### Download repository

```
git clone https://github.com/fanzhuyifan/Nanopore-Compression.git
```

### Install
The instructions below will create the executable in the build directory inside Nanopore-Compression.

On Linux with cmake installed and version at least 3.10 (check using cmake --version):
```
cd Nanopore-Compression
mkdir build
cd build
cmake ..
make
```
## Usage
Run the Nanopore-Compression executable with the options below:
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
Note that the compressed files are tar archives consisting of the different compressed streams, although we recommend using the .nanopore-compression extension as in the examples shown below.

## Example Usage of Nanopore Compression
This section contains several examples for compression and decompression with various modes and options. The compressed file uses the .nanopore-compression extension as a convention.

For compressing file.fastq losslessly using default 16 threads (Lossless).
```
./testCompressor -i file.fastq -o file.nanopore-compression
```
Using 20 threads (Lossless).
```
./testCompressor -i file.fastq -o file.nanopore-compression -t 20
```
For compressing file.fastq.gz (gzipped fastq files) losslessly using default 16 threads (Lossless).
```
./testCompressor -i file.fastq -o file.nanopore-compression -g 
```
Compressing with kmer size for minhash 23, number of hash functions for minhash 54, and the overlap sketch threhold for minhash 6.
```
./testCompressor -i file.fastq -k 23 -n 54 --overlap-sketch-thr 6 -o file.nanopore-compression  
```
Compressing with kmer size of minimap2 25 and the window size for the minimap2 100.
```
./testCompressor -i file.fastq --minimap-k 25 --minimap-w 100 -o file.nanopore-compression  
```
Compressing with the max number of partial chains during chaining for minimap2 to be 800.
```
./testCompressor -i file.fastq --max-chain-iter 800 -o file.nanopore-compression  
```
Compressing with the max number of edges allowed in a consensus graph to be 4000000.
```
./testCompressor -i file.fastq --edge-thr 4000000 -o file.nanopore-compression  
```
Decompressing to file.reads.
```
./testDecompressor file.nanopore-compression file.reads
```
