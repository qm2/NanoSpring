#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v ./../../../minimap2-2.18/minimap2 -x map-ont --secondary=no --cs $1".fasta" $1".fastq" -t 20 > $1".paf"  |& tee -a $2

stdbuf -oL /usr/bin/time -v renano -t 20 -s $1".fasta" $1".paf" $1".fastq" $1".renano" |& tee -a $2