#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"nt50_4.NanoSpring" -n 50 --overlap-sketch-thr 4 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"nt50_5.NanoSpring" -n 50 --overlap-sketch-thr 5 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"nt60_6.NanoSpring" -n 60 --overlap-sketch-thr 6 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"nt70_7.NanoSpring" -n 70 --overlap-sketch-thr 7 |& tee -a $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"nt50_6.NanoSpring" -n 50 --overlap-sketch-thr 6 |& tee -a $2