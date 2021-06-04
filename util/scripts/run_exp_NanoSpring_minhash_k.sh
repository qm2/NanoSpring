#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"k30.NanoSpring" -k 30 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"k25.NanoSpring" -k 25 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"k23.NanoSpring" -k 23 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"k20.NanoSpring" -k 20 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"k18.NanoSpring" -k 18 |& tee -a $2