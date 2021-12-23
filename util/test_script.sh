#!/bin/bash
# exit when any command fails
set -e

echo "Running test"

./NanoSpring -c -i ../util/test_file.fastq.gz -o tmp.nanospring
./NanoSpring -d -i tmp.nanospring -o tmp.reads
cmp <(zcat ../util/test_file.fastq.gz | sed -n '2~4p') tmp.reads 
rm tmp.nanospring
rm tmp.reads

./NanoSpring -c -i ../util/test_file.fastq.gz -o tmp.nanospring --low-memory-compression
./NanoSpring -d -i tmp.nanospring -o tmp.reads
cmp <(zcat ../util/test_file.fastq.gz | sed -n '2~4p') tmp.reads 
rm tmp.nanospring
rm tmp.reads

echo "Test was successful"
