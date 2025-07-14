import re
import sys

def parse_alignment_file(filename):
    """Extracts the total number of reads from an alignment summary file."""
    total_reads = 0

    with open(filename, 'r') as f:
        for line in f:
            if match := re.search(r'(\d+) reads;', line):
                total_reads = int(match.group(1))
                break  # Stop after finding the first match

    return total_reads

def parse_patterns_file(filename):
    """Extracts the number of aligned reads from a patterns file."""
    aligned_reads = 0

    with open(filename, 'r') as f:
        for line in f:
            if match := re.search(r'(\d+)\s+patterns loaded from file', line):
                aligned_reads = int(match.group(1))
                break  # Stop after finding the first match

    return aligned_reads

def process_files(file_pairs, output_file):
    """Processes multiple pairs of files and writes the summed results to an output file."""
    total_reads_sum = 0
    uniquely_aligned_sum_total = 0

    for i in range(0, len(file_pairs), 2):
        alignment_file = file_pairs[i]
        patterns_file = file_pairs[i + 1]

        total_reads = parse_alignment_file(alignment_file)
        uniquely_aligned_reads = parse_patterns_file(patterns_file)

        total_reads_sum += total_reads
        uniquely_aligned_sum_total += uniquely_aligned_reads

    # Calculate percentage
    percentage = (uniquely_aligned_sum_total / total_reads_sum) * 100 if total_reads_sum > 0 else 0

    # Write results to output file
    with open(output_file, 'w') as out:
        out.write(f"{uniquely_aligned_sum_total} aligned to exclusion index out of {total_reads_sum} reads.\n")
        out.write(f"Percent excluded: {percentage:.2f}%\n")

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) % 2 == 1:
        print("Usage: python merge_exclusion_stats.py <output_file> [input file pairs] ...")
        sys.exit(1)

    output_filename = sys.argv[1]
    input_file_pairs = sys.argv[2:]

    process_files(input_file_pairs, output_filename)

