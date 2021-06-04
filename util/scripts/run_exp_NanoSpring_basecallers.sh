#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1"_fast.fastq" -o $1"_fast.NanoSpring" |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1"_hac.fastq" -o $1"_hac.NanoSpring" |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1"_sup.fastq" -o $1"_sup.NanoSpring" |& tee -a $2