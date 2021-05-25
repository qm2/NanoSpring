#!/bin/bash
echo "The input fastq file path is:" $1
echo "The log file path is:" $2

stdbuf -oL /usr/bin/time -v ./../../build/NanoSpring -c -i $1 -o $1".NanoSpring" |& tee -a $2
stdbuf -oL /usr/bin/time -v ./../../build/NanoSpring -d -i $1".NanoSpring" -o $1".NanoSpring.reads" |& tee -a $2
stdbuf -oL cmp $1".NanoSpring.reads" <(sed -n '2~4p' $1) |& tee -a $2

stdbuf -oL /usr/bin/time -v ./../../../EnanoFASTQ/enano/enano -t 20 $1 $1".enano" |& tee -a $2
stdbuf -oL /usr/bin/time -v ./../../../EnanoFASTQ/enano/enano -d -t 20 $1".enano" enano_decompressed.fastq |& tee -a $2
stdbuf -oL cmp $1 enano_decompressed.fastq |& tee -a $2
rm enano_decompressed.fastq

stdbuf -oL /usr/bin/time -v pigz -p 20 $1".NanoSpring.reads" |& tee -a $2
stdbuf -oL ls -l $1".NanoSpring.reads.gz"
stdbuf -oL /usr/bin/time -v pigz -d -k -p 20 $1".NanoSpring.reads.gz" |& tee -a $2

rm $1".NanoSpring.reads"

