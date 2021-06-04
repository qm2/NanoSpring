#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v /raid/qingxi/NanoSpring/build/NanoSpring -c -i $1 -o $1"_edge_unlimited.NanoSpring" |& tee -a $2