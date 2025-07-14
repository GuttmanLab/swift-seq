import pysam
import sys

# Prepares a "hybrid genotyped" VCF where the two alleles represent two different species
# Essentially, the genotype (GT) field is reformatted to x/y where x is the allele of species 1 and y is the allele of species 2
# The input VCF file should contain both species across two samples columns

# Arguments: <input_vcf_file> <sample_1> [<sample_2>] <output_vcf_file>

def main():

    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python " + sys.argv[0] + " <input_vcf_file> <sample_1> [<sample_2>] <output_vcf_file>")
        sys.exit(1)

    s1 = sys.argv[2]
    s2 = sys.argv[3] if len(sys.argv) == 5 else ""

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[4] if s2 else sys.argv[3]

    # Handle cases where sample_1 or sample_2 is "" or "REF"
    if s1 == "" or s1.upper() == "REF":
        s1 = s2
        s2 = ""
    if s2 == "" or s2.upper() == "REF":
        s2 = None

    if s1 == "" or s1.upper() == "REF":
        print("Error: At least one of sample_1 or sample_2 must be specified and not be REF or empty.")
        sys.exit(1)

    vcf_in = pysam.VariantFile(input_vcf, "r")
    header = vcf_in.header
    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)

    # Prepare some logs
    log_s1_not_biallelic = 0
    log_s2_not_biallelic = 0
    log_s1_not_homo = 0
    log_s2_not_homo = 0
    log_s1_s2_identical = 0
    log_n = 0
    log_n_read = 0

    for rec in vcf_in:
        # Process each record
        s1_gt = rec.samples[s1]['GT']  # Genotype of sample 1
        s2_gt = rec.samples[s2]['GT'] if s2 else None  # Genotype of sample 2
        log_n_read += 1

        # Create new genotype field based on these indices
        # Check whether there are exactly 2 alleles (if not, discard)
        if len(s1_gt) != 2:
            log_s1_not_biallelic += 1
            continue
        if s2 and len(s2_gt) != 2:
            log_s2_not_biallelic += 1
            continue
        # Check whether the genotypes are homozygous (if not, discard)
        if s1_gt[0] != s1_gt[1]:
            log_s1_not_homo += 1
            continue
        if s2 and s2_gt[0] != s2_gt[1]:
            log_s2_not_homo += 1
            continue
        # Make the "new" genotype
        s1_new_gt = s1_gt[0]
        s2_new_gt = s2_gt[0] if s2 else s1_gt[0]  # Use s1_gt[0] for single species mode
        if not s2:
            s1_new_gt = "0" # Use 0 for single species mode
        # Replace . with 0
        if s1_new_gt == "." or s1_new_gt == None:
            s1_new_gt = "0"
        if s2_new_gt == "." or s2_new_gt == None:
            s2_new_gt = "0"

        # Identify homozygous genotypes
        if s1_new_gt == s2_new_gt:
            log_s1_s2_identical += 1

        new_gt = (s1_new_gt, s2_new_gt)
        log_n += 1

        # Set new genotype to all samples
        for sample in rec.samples:
            rec.samples[sample]['GT'] = new_gt

        # Write the modified record
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()

    # Print out the log
    print(f'Sample 1: {s1}, Sample 2: {s2 if s2 else "None (single species mode)"}')
    print(f'Read in {log_n_read} records from input VCF file')
    print(f'Discarded {log_s1_not_biallelic} records because not exactly two alleles in sample {s1}')
    if s2:
        print(f'Discarded {log_s2_not_biallelic} records because not exactly two alleles in sample {s2}')
    print(f'Discarded {log_s1_not_homo} records because not homozygous in sample {s1}')
    if s2:
        print(f'Discarded {log_s2_not_homo} records because not homozygous in sample {s2}')
    print(f'Note: Output VCF file has {log_s1_s2_identical} homozygous genotypes')
    print(f'Wrote out {log_n} records into output VCF file')

if __name__ == "__main__":
    main()


