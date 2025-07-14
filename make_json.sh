#!/bin/bash

set -e

declare -A sample_R1
declare -A sample_R2

for dir in "$@"; do
  for fq in "$dir"/*_R[12].fastq.gz; do
    [[ -e "$fq" ]] || continue  # skip if glob doesn't match anything
    fq_base=$(basename "$fq")
    fq_path="$dir/$fq_base"
    # Extract sample name before first _R1 or _R2
    sample=$(echo "$fq_base" | sed -E 's/_R[12].fastq.gz$//' | sed -E 's/[_-]?[0-9]+$//' | sed -E 's/[-_]+$//')
    [[ -z "$sample" ]] && continue
    if [[ "$fq_base" == *_R1.fastq.gz ]]; then
      sample_R1["$sample"]+="$fq_path,"
    elif [[ "$fq_base" == *_R2.fastq.gz ]]; then
      sample_R2["$sample"]+="$fq_path,"
    fi
  done
done

# Write JSON output
json_file="./example_samples.json"
echo "{" > "$json_file"
samples=("${!sample_R1[@]}" "${!sample_R2[@]}")
samples=($(printf "%s\n" "${samples[@]}" | sort -u))
last_index=$((${#samples[@]} - 1))
for i in "${!samples[@]}"; do
  sample="${samples[$i]}"
  r1_files="[\"$(echo "${sample_R1[$sample]}" | sed 's/,$//' | sed 's/,/\", \"/g')\"]"
  r2_files="[\"$(echo "${sample_R2[$sample]}" | sed 's/,$//' | sed 's/,/\", \"/g')\"]"
  echo "  \"$sample\": {" >> "$json_file"
  echo "    \"R1\": $r1_files," >> "$json_file"
  echo "    \"R2\": $r2_files" >> "$json_file"
  if [[ "$i" -eq "$last_index" ]]; then
    echo "  }" >> "$json_file"
  else
    echo "  }," >> "$json_file"
  fi
done
echo "}" >> "$json_file"

