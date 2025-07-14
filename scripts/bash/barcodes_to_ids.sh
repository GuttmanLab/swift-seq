#!/bin/bash

# Converts the synthetically generated barcodes (e.g. AAAAAAAAAAAAAAGG) to a barcode ID, e.g. [NYBot13_Stg][Odd2Bo28][RTBC7], based on a mapping file 

# Check if exactly two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <mapping_file> <barcodes_file>"
    exit 1
fi

# Check if the first argument is a file
if [ ! -f "$1" ]; then
    echo "Error: '$1' is not a file or does not exist."
    exit 1
fi

# Check if the second argument is a file
if [ ! -f "$2" ]; then
    echo "Error: '$2' is not a file or does not exist."
    exit 1
fi

mapping_file="$1"
barcodes_file="$2"

awk 'NR==FNR {arr[$1] = $2; next} {if ($1 in arr) {split(arr[$1], a, ","); for (i in a) printf "[%s]%s", a[i], (i==length(a) ? "\n" : "")} }' "$mapping_file" "$barcodes_file"


