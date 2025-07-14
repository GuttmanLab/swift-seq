#!/usr/bin/env python3

import argparse
import re
import os
from collections import defaultdict

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Concatenate ligation efficiency files by summing counts and recalculating percentages."
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="Path to the final concatenated output file."
    )
    parser.add_argument(
        'input_files',
        nargs='+',
        help="One or more input ligation efficiency files to concatenate."
    )
    return parser.parse_args()

def parse_file(file_path, reads_with_barcodes, barcodes_in_position):
    """
    Parses a single input file and updates the counts in the provided dictionaries.

    Args:
        file_path (str): Path to the input file.
        reads_with_barcodes (defaultdict): Dictionary to accumulate reads with X barcodes.
        barcodes_in_position (defaultdict): Dictionary to accumulate barcodes found in position Y.
    """
    # Define regex patterns
    reads_pattern = re.compile(r'^(\d+)\s+\(([\d.]+)%\)\s+reads found with\s+(\d+)\s+barcode[s]?\.?$')
    barcodes_pattern = re.compile(r'^(\d+)\s+\(([\d.]+)%\)\s+barcodes found in position\s+(\d+)\.?$')

    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            # Try matching reads with barcodes
            reads_match = reads_pattern.match(line)
            if reads_match:
                count, percent, x_barcodes = reads_match.groups()
                x_barcodes = int(x_barcodes)
                count = int(count)
                reads_with_barcodes[x_barcodes] += count
                continue

            # Try matching barcodes in position
            barcodes_match = barcodes_pattern.match(line)
            if barcodes_match:
                count, percent, position = barcodes_match.groups()
                position = int(position)
                count = int(count)
                barcodes_in_position[position] += count
                continue

            # If line doesn't match any pattern, issue a warning
            print(f"Warning: Line {line_num} in file '{file_path}' does not match expected patterns:\n  {line}")

def concatenate_files(output_path, reads_with_barcodes, barcodes_in_position):
    """
    Writes the aggregated counts and recalculated percentages to the output file.

    Args:
        output_path (str): Path to the output file.
        reads_with_barcodes (dict): Aggregated counts for reads with X barcodes.
        barcodes_in_position (dict): Aggregated counts for barcodes found in position Y.
    """
    # Calculate total counts
    total_reads = sum(reads_with_barcodes.values())
    total_barcodes = sum(barcodes_in_position.values())

    # Sort the keys for consistent ordering
    sorted_barcodes = sorted(reads_with_barcodes.keys())
    sorted_positions = sorted(barcodes_in_position.keys())

    with open(output_path, 'w') as out:
        # Write Reads Found with X Barcodes
        for x in sorted_barcodes:
            count = reads_with_barcodes[x]
            percentage = (count / total_reads) * 100 if total_reads > 0 else 0
            out.write(f"{count} ({percentage:.1f}%) reads found with {x} barcode{'s' if x != 1 else ''}.\n")

        out.write("\n")  # Add a blank line between sections

        # Write Barcodes Found in Position Y
        for y in sorted_positions:
            count = barcodes_in_position[y]
            percentage = (count / total_reads) * 100 if total_reads > 0 else 0
            out.write(f"{count} ({percentage:.1f}%) barcodes found in position {y}.\n")

def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Initialize dictionaries to hold aggregated counts
    reads_with_barcodes = defaultdict(int)    # key: X barcodes, value: total count
    barcodes_in_position = defaultdict(int)   # key: Y position, value: total count

    # Process each input file
    for file_path in args.input_files:
        if not os.path.isfile(file_path):
            print(f"Error: Input file '{file_path}' does not exist or is not a file.")
            continue
        parse_file(file_path, reads_with_barcodes, barcodes_in_position)

    # Generate the concatenated output
    concatenate_files(args.output, reads_with_barcodes, barcodes_in_position)

    print(f"Concatenation complete. Output written to '{args.output}'.")

if __name__ == "__main__":
    main()

