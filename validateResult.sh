#!/bin/bash
diff -s <(sed '1~2d' "$1") "$2"
if [[ $? -eq 0 ]]; then
    echo "Succeeded"
else
    echo "Failed"
fi
