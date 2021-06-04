#!/bin/bash
'''
based on https://github.com/hcdenbakker/N50.sh/blob/master/N50.sh
Obtain statistics for fastq file
Usage
cat file.fastq | sed -n '2~4p' | awk '{ print length }' | ./N50.sh
zcat file.fastq.gz | sed -n '2~4p' | awk '{ print length }' | ./N50.sh
'''

mkdir temp;
#get read lengths, ordered from large to small
cat - | sort -gr > temp/read_lengths.txt;

#number of reads
Y=$(cat temp/read_lengths.txt | wc -l);

#sum of read_lengths
X=$(paste -sd+ temp/read_lengths.txt | bc);

# cumulative read lengths

awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' temp/read_lengths.txt > temp/read_lengths_cum.txt;

# get cumulative read contributions (%) to the entire file

awk -v var=$X 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) $i/=var;}1' temp/read_lengths_cum.txt > temp/cum_perc.txt; 

# join results

paste temp/read_lengths.txt temp/cum_perc.txt > temp/matrix.txt;

# get N50, N90, largest scaffold/read

N50=$(awk '$2 >= 0.50' temp/matrix.txt |head -1| awk '{ print $1}');
N90=$(awk '$2 >= 0.90' temp/matrix.txt |head -1| awk '{ print $1}');
large_read=$(head -1 temp/read_lengths.txt);
rm -r temp;

echo "number of reads:$Y"
echo "total size:$X"
echo "largest read:$large_read"
echo "N50:$N50";
echo "N90:$N90";
