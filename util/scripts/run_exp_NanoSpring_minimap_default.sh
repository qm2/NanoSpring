#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2
# gzip -d -k $1".gz"

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"minimap_1.NanoSpring" --minimap-k 15 --minimap-w 10 --max-chain-iter 5000  |& tee -a $2

# stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"minimap_2.NanoSpring" --minimap-k 20 --minimap-w 50 --max-chain-iter 5000  |& tee -a $2