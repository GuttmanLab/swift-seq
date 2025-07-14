# Usage: python generate_barcodes.py <barcode_length> <number_of_barcodes> <pre_existing_barcodes_file [optional]>
# This script generates barcodes, greedily maximizing the hamming distances between each barcode
# For example, `python generate_barcodes.py 5 96 pre_existing_barcodes.txt` would output 96 5-bp barcodes with pre_existing_barcodes.txt being the barcodes that we already have
# Author: Delaney K. Sullivan

from itertools import product
import sys
import random


def hamming_distance(seq1, seq2):
	"""Calculate the Hamming distance between two sequences."""
	return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def generate_all_sequences(length=5):
	"""Generate all possible sequences of given length."""
	nucleotides = ['A', 'T', 'C', 'G']
	return [''.join(p) for p in product(nucleotides, repeat=length)]


def filter_sequences(initial_sequences, all_sequences, min_distance=2, target_count=96):
	"""Filter sequences to ensure minimum Hamming distance, including specified sequences."""
	filtered_sequences = initial_sequences[:]
	for seq in all_sequences:
		if seq in initial_sequences:  # Skip if sequence is already included
			continue
		if all(hamming_distance(seq, existing_seq) >= min_distance for existing_seq in filtered_sequences):
			filtered_sequences.append(seq)
			if len(filtered_sequences) == target_count:
				break
	return filtered_sequences



if len(sys.argv) < 3:
	print("First argument: Barcode length")
	print("Second argument: Number of barcodes")
	print("Third argument [optional]: File with pre-existing barcodes")
	sys.exit(1)

lines = []
if len(sys.argv) >= 4:
	with open(sys.argv[3]) as file:
	    	lines = [line.rstrip() for line in file]
specified_sequences = lines

if len(specified_sequences) >= int(sys.argv[2]):
	print("Error: Number of pre-existing barcodes cannot exceed number of barcodes requested")
	sys.exit(1)


all_sequences = generate_all_sequences(int(sys.argv[1]))
random.shuffle(all_sequences)

dist = int(sys.argv[1])
num = int(sys.argv[2])

while dist >= 0:
	filtered_sequences = filter_sequences(specified_sequences, all_sequences, dist, num)
	if len(filtered_sequences) == num:
		break
	else:
		dist = dist - 1
		specified_sequences = filtered_sequences.copy()
	
for seq in filtered_sequences:
	print(seq)

