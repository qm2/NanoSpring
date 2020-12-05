#!/bin/bash
# Usage: ./validateResult.sh originalFile decompressedFile 
# originalFile must end in .reads or .fastq

ext="${1##*.}"
if [ "$ext" = "reads" ]; then
    diff -s <(sed '1~2d' "$1") "$2"
elif [ "$ext" = "fastq" ]; then
    diff -s <(sed '1~4d;3~4d;4~4d;' "$1") "$2"
else
    echo "Unknown extension: $1"
    exit
fi

if [[ $? -eq 0 ]]; then
    echo "Succeeded"
else
    echo "Failed"
fi
