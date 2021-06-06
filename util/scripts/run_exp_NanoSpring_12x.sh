#!/bin/bash
file1="/raid/qingxi/Nanopore-Compression-old/data/ERR5455028_13x.fastq"
file2="/raid/qingxi/Nanopore-Compression-old/data/rel7_sampled_13x.fastq"
logname="/raid/qingxi/NanoSpring/logs/13x_coverage_experiment.log"

stdbuf -oL /usr/bin/time -v ./../../build/NanoSpring -c -i $file1 -o $file1".NanoSpring" |& tee -a $logname

stdbuf -oL /usr/bin/time -v ./../../build/NanoSpring -c -i $file2 -o $file2".NanoSpring" |& tee -a $logname