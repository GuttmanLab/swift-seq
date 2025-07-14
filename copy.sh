#!/bin/bash

set -e

# Check at least one argument is given
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <target_dir> [input_dir1 input_dir2 ...]"
  exit 1
fi

# First argument: destination parent directory
base_dest="$1"
shift
target_dir="$base_dest/swiftseq_pipeline"

# Ensure target_dir does not already exist
if [[ -e "$target_dir" ]]; then
  echo "Error: Target directory '$target_dir' already exists."
  exit 1
fi

# Create the target directory
mkdir -p "$target_dir"

# Copy contents of current directory into target_dir
cp -r ./* "$target_dir"

# If there are additional arguments, verify and process them
if [[ $# -gt 0 ]]; then
  abs_dirs=()
  for d in "$@"; do
    if [[ ! -d "$d" ]]; then
      echo "Error: Directory '$d' does not exist."
      exit 1
    fi
    fq_count=$(find "$d" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)
    if [[ "$fq_count" -eq 0 ]]; then
      echo "Error: Directory '$d' does not contain any FASTQ files."
      exit 1
    fi
    abs_dirs+=("$(cd "$d" && pwd)")
  done

  # Call make_json.sh with absolute directories
  "$target_dir/make_json.sh" "${abs_dirs[@]}"
fi

# Change into the target directory
cd "$target_dir"

ulimit -n 8192

