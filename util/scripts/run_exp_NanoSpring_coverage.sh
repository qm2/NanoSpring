#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring-test/build/NanoSpring -c -i $1"_25x.fastq" -o $1"_25x.fastq.NanoSpring" |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring-test/build/NanoSpring -c -i $1"_50x.fastq" -o $1"_50x.fastq.NanoSpring" |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring-test/build/NanoSpring -c -i $1"_100x.fastq" -o $1"_100x.fastq.NanoSpring" |& tee -a $2
