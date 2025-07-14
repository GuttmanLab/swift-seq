'''
Author: Delaney K. Sullivan
Aim: A Snakemake workflow to process SWIFT-seq data
'''

import json
import re
import os 
import sys
import datetime
import gzip
from os.path import basename, splitext
from pathlib import Path

wildcard_constraints:
    chunk="(?:\.part_.*|[^.].*)?" # Permits us to have empty string as a wildcard; disallows leading periods w/ exception of .part_

################################################################################
#Load config.yaml file and other general settings
################################################################################

#Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')


try:
    config_path = config["config_path"]
except:
    config_path = 'config.yaml'

configfile: config_path

try:
    num_chunks = int(config["num_chunks"])
    if num_chunks <= 0:
        num_chunks = 1
    print("Number of chunks specified:", num_chunks, file=sys.stderr)
except:
    num_chunks = 1
    print("No valid 'num_chunks' in config file. Defaulting to 1 (no chunking).", file=sys.stderr)

try:
    keep_fastq_chunks = config["keep_fastq_chunks"]
    if num_chunks > 1 and keep_fastq_chunks:
        print("Will keep chunked FASTQ files", file=sys.stderr)
    else:
        keep_fastq_chunks = False
except:
    keep_fastq_chunks = False

try:
    chunk_star_alignment = config["chunk_star_alignment"]
    if num_chunks > 1 and chunk_star_alignment:
        print("Using chunks when performing STAR alignment", file=sys.stderr)
    else:
        chunk_star_alignment = False
except:
    chunk_star_alignment = False

try:
    thread_factor = int(config["thread_factor"])
    print("Thread factor specified:", num_chunks, file=sys.stderr)
except:
    thread_factor = 1
    print("No valid 'thread_factor' in config file. Defaulting to 1.", file=sys.stderr)

try:
    email = config['email']
except:
    print("Won't send email on error")
    email = None
try:
    samples = os.path.abspath(config["samples"])
    print("Using samples file: ", samples, file=sys.stderr)
except:
    samples = "./samples.json"
    print("Defaulting to samples JSON file: ", samples, file=sys.stderr)
try:
    R1_adapter = config["R1_adapter"]
    R2_adapter = config["R2_adapter"]
    print("Using R1 adapter: ", R1_adapter, file=sys.stderr)
    print("Using R2 adapter: ", R2_adapter, file=sys.stderr)
except:
    R1_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    R2_adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    print("Defaulting to R1 adapter: ", R1_adapter, file=sys.stderr)
    print("Defaulting to R2 adapter: ", R2_adapter, file=sys.stderr)
try:
    splitcode_bID = os.path.abspath(config["bID"])
    print("Using splitcode config file: ", splitcode_bID, file=sys.stderr)
except:
    splitcode_bID = "config.txt"
    print("Defaulting to splitcode file: ", splitcode_bID, file=sys.stderr)
umi_extraction_pattern = "0:0<R1>0:-1,<R1[10]>{adapter},{adapter}<R1[1-65]>,2:0<R2[1-65]>"
try:
    if config["umi_extraction_pattern"].strip():
        umi_extraction_pattern = config["umi_extraction_pattern"].strip()
    print("Using UMI extraction pattern: ", umi_extraction_pattern, file=sys.stderr)
except:
    print("Defaulting to UMI extraction pattern: ", umi_extraction_pattern, file=sys.stderr)
try:
    tag_prefix_order=config["tag_prefix_order"]
    print("Using tag prefix order: ", tag_prefix_order, file=sys.stderr)
except:
    tag_prefix_order=""
    print("Skipping calculating ligation efficiency", file=sys.stderr)
try:
    temp_dir = os.path.abspath(config["temp_dir"])
    print("Using temporary directory: ", temp_dir, file=sys.stderr)
except:
    temp_dir = "tmp/"
    print("Defaulting to temporary directory: ", temp_dir, file=sys.stderr)
try:
    subset_n = int(config["subset_n"])
    if subset_n > 0:
        print("SUBSETTING: ", str(subset_n), " reads", file=sys.stderr)
    else:
        subset_n = 0
except:
    subset_n = 0
try:
    fastqs_dir = config["fastqs"]
    if not fastqs_dir or fastqs_dir == "./":
        fastqs_dir=""
        print("Location of FASTQS file is the current directory", file=sys.stderr)
    else:
        print("Location of FASTQS file: ", temp_dir, file=sys.stderr)
except:
    fastqs_dir = ""
    print("Default to current directory for location of FASTQS file: ", file=sys.stderr)
try:
    align_id = config["align_id"]
    print("Alignment ID: ", align_id, file=sys.stderr)
except:
    align_id="default"
    print("Using default alignment ID", file=sys.stderr)
try:
    genome_fasta = config["genome_fasta"]
    genome_gtf = config["genome_gtf"]
except:
    genome_fasta = ""
    genome_gtf = ""
try:
    use_features = config["use_features"]
except:
    use_Features = False
try:
    features_fasta_file = ""
    if use_features:
        features_fasta_file = config["features_fasta_file"]
        if features_fasta_file:
            print("Using features: ", features_fasta_file, file=sys.stderr)
        else:
            use_features = False
except:
    use_features = False
try:
    features_kmer_length = 0
    if use_features:
        features_kmer_length = int(config["features_kmer_length"])
        if not features_kmer_length or features_kmer_length < 3:
            features_kmer_length = 31
            print("Invalid features k-mer length supplied, defaulting to ", features_kmer_length, file=sys.stderr)
        else:
            print("Features k-mer length set to ", features_kmer_length, file=sys.stderr)
except:
    if use_features:
        features_kmer_length = 31
        print("No features k-mer length supplied, defaulting to ", features_kmer_length, file=sys.stderr)
try:
    features_start = 0
    if use_features:
        features_start = int(config["features_start"])
        if not features_start or features_start <= 0:
            features_start = 0
        else:
            print("Features start position set to ", features_start, file=sys.stderr)
except:
    features_start = 0
try:
    features_length = -1
    if use_features:
        features_length = int(config["features_length"])
        if not features_length or features_length <= 0:
            features_length = -1
        else:
            print("Features sequence length set to ", features_length, file=sys.stderr)
except:
    features_length = -1

try:
    use_vcf_wasp_workflow = config["use_vcf_wasp_workflow"]
    if use_vcf_wasp_workflow:
        print("Variant-specific mapping (WASP) enabled", file=sys.stderr)
except:
    use_vcf_wasp_workflow = False
try:
    vcf_file = config["vcf_file"]
    print("VCF file: ", vcf_file, file=sys.stderr)
    if not use_vcf_wasp_workflow:
        print("Not using VCF file", file=sys.stderr)
        vcf_file = ""
except:
    vcf_file = ""
try:
    vcf_s1 = config["vcf_s1"]
    vcf_s2 = config["vcf_s2"]
except:
    vcf_s1 = ""
    vcf_s2 = ""
vcf_s1_ = "REF"
if vcf_file and use_vcf_wasp_workflow:
    if vcf_s1:
        vcf_s1_ = vcf_s1
    print("VCF sample 1: ", vcf_s1_, file=sys.stderr)
    if not vcf_s2:
        print("Must supply VCF sample 2 for variant-specific mapping workflow", file=sys.stderr)
        sys.exit()
    else:
        print("VCF sample 2: ", vcf_s2, file=sys.stderr)
try:
    split_allele_analysis = config["split_allele_analysis"]
    if split_allele_analysis and vcf_file and use_vcf_wasp_workflow:
        print("Turning on split allele analysis", file=sys.stderr)
    else:
        split_allele_analysis = False
except:
    split_allele_analysis = False

try:
    use_exclusion_filter = config["use_exclusion_filter"]
    exclusion_fasta_file = ""
    if use_exclusion_filter:
        exclusion_fasta_file = config["exclusion_fasta_file"]
        if exclusion_fasta_file:
            print("Enabling exclusion-filtering for reads that align to file: ", exclusion_fasta_file, file=sys.stderr)
        else:
            use_exclusion_filter = False
except:
    use_exclusion_filter = False
    exclusion_fasta_file = ""


# STAR settings:
try:
    use_star = False
    if config["use_star"]:
        print("STAR aligner enabled", file=sys.stderr)
        use_star=True
except:
    use_star = False
try:
    star_index = config["star_index"]
except:
    star_index = None
if use_star:
    if star_index:
        star_index = os.path.abspath(star_index)
        print("Using provided STAR index: ", star_index, file=sys.stderr)
    elif genome_fasta and genome_gtf:
        genome_fasta = os.path.abspath(genome_fasta)
        genome_gtf = os.path.abspath(genome_gtf)
        print("STAR index generation from FASTA and GTF: ", genome_fasta, " ", genome_gtf, file=sys.stderr)
    else:
        print("Disabling STAR because index not supplied, and FASTA+GTF also not supplied", file=sys.stderr)
        use_star = False
try:
    star_additional_params = ""
    if config["use_star"] and config["star_additional_params"]:
        star_additional_params = config["star_additional_params"]
        print("STAR aligner set to run with additional parameters: ", star_additional_params, file=sys.stderr)
except:
    star_additional_params = ""


# kallisto settings:
try:
    use_kallisto = False
    if config["use_kallisto"]:
        print("kallisto enabled", file=sys.stderr)
        use_kallisto=True
except:
    use_kallisto = False
try:
    kallisto_index = config["kallisto_index"]
    kallisto_index_t2g = config["kallisto_index_t2g"]
    kallisto_index_c1 = config["kallisto_index_c1"]
    kallisto_index_c2 = config["kallisto_index_c2"]
    if not (kallisto_index and kallisto_index_t2g and kallisto_index_c1 and kallisto_index_c2):
        kallisto_index = None
        kallisto_index_t2g = None
        kallisto_index_c1 = None
        kallisto_index_c2 = None
    else:
        kallisto_index = os.path.abspath(kallisto_index)
        kallisto_index_t2g = os.path.abspath(kallisto_index_t2g)
        kallisto_index_c1 = os.path.abspath(kallisto_index_c1)
        kallisto_index_c2 = os.path.abspath(kallisto_index_c2)
except:
    kallisto_index = None
    kallisto_index_t2g = None
    kallisto_index_c1 = None
    kallisto_index_c2 = None

if use_kallisto:
    if kallisto_index:
        print("Using provided kallisto index: ", kallisto_index, file=sys.stderr)
        print("Using provided kallisto transcripts-to-gene mapping file: ", kallisto_index_t2g, file=sys.stderr)
    elif genome_fasta and genome_gtf:
        genome_fasta = os.path.abspath(genome_fasta)
        genome_gtf = os.path.abspath(genome_gtf)
        print("kallisto index generation from FASTA and GTF: ", genome_fasta, " ", genome_gtf, file=sys.stderr)
    else:
        print("Disabling kallisto because index not supplied, and FASTA+GTF also not supplied", file=sys.stderr)
        use_kallisto = False
try:
    kallisto_dlist = config["kallisto_dlist"]
    if use_kallisto and not kallisto_index:
        if kallisto_dlist and kallisto_dlist.upper() == "NONE":
            kallisto_dlist = "None"
            print("kallisto index d-list disabled", file=sys.stderr)
        elif kallisto_dlist:
            kallisto_dlist = ','.join([os.path.abspath(x) for x in kallisto_dlist.split(",")]) # Allow multiple files
            print("kallisto index d-list set to: ", kallisto_dlist, file=sys.stderr)
        else:
            kallisto_dlist = genome_fasta
            print("kallisto index d-list not provided; defaulting to genome FASTA provided", file=sys.stderr)
except:
    kallisto_dlist = None
    if use_kallisto and not kallisto_index:
        kallisto_dlist = genome_fasta
        print("kallisto index d-list not provided; defaulting to genome FASTA provided", file=sys.stderr)



# Other settings:
try:
    select_tags = ""
    if config["select_tags"]:
        select_tags = config["select_tags"]
        print("Only including barcodes with the following tags in the final AnnData: ", select_tags, file=sys.stderr)
except:
    select_tags = ""
try:
    merge_anndata_files = False
    if config["merge_anndata"]:
        merge_anndata_files = True
        print("Outputting a merged anndata object file is enabled", file=sys.stderr)
except:
    merge_anndata_files = False
try:
    out_dir = config["output_dir"]
    print("All data will be written to: ", out_dir, file=sys.stderr)
except:
    out_dir = os.getcwd()
    print("Defaulting to working directory as output directory", file=sys.stderr)
try:
    conda_env = config["conda_env"]
except:
    print("No conda environment specified. Defaulting to envs/swiftseq.yaml", file=sys.stderr)
    conda_env = "envs/swiftseq.yaml"
if conda_env.lower().endswith(".yaml") or conda_env.lower().endswith(".yml"):
    print("Will create new conda environment from", conda_env, file=sys.stderr)
else:
    print("Using existing conda environment:", conda_env, file=sys.stderr)

adata_filename="adata.h5ad"

##############################################################################
# Location of scripts
##############################################################################

try:
    DIR_SCRIPTS = config["scripts_dir"]
except:
    print("Scripts directory not specificed in config.yaml", file=sys.stderr)
    sys.exit()  # no default, exit

split_fastq = os.path.join(DIR_SCRIPTS, "bash/split_fastq.sh")
ligeff = os.path.join(DIR_SCRIPTS, "bash/ligeff.sh")
merge_ligeff = os.path.join(DIR_SCRIPTS, "python/merge_ligeff.py")
merge_exclusion_stats = os.path.join(DIR_SCRIPTS, "python/merge_exclusion_stats.py")
barcodes_to_ids  = os.path.join(DIR_SCRIPTS, "bash/barcodes_to_ids.sh")
prep_anndata = os.path.join(DIR_SCRIPTS, "python/prep_anndata.py")
merge_anndata = os.path.join(DIR_SCRIPTS, "python/merge_anndata.py")
prep_vcf = os.path.join(DIR_SCRIPTS, "python/prep_vcf.py")
decode_barcode = os.path.join(DIR_SCRIPTS, "python/decode_barcode.py")
split_bam_by_snp = os.path.join(DIR_SCRIPTS, "java/SplitBAMBySNP.jar")
html_report_gen = os.path.join(DIR_SCRIPTS, "python/html/report.py")

################################################################################
#make output directories (aren't created automatically on cluster)
################################################################################

DIR_WORKUP = os.path.join(out_dir, "workup")
DIR_LOGS = os.path.join(DIR_WORKUP, "logs")

DIR_LOGS_CLUSTER = os.path.join(DIR_LOGS, "cluster")
os.makedirs(DIR_LOGS_CLUSTER, exist_ok=True)
out_created = os.path.exists(DIR_LOGS_CLUSTER)
print("Output logs path created:", out_created, file=sys.stderr)

DIR_TRIMMED_CLUSTER = os.path.join(DIR_WORKUP, "trimmed", "cluster")
os.makedirs(DIR_TRIMMED_CLUSTER, exist_ok=True)
out_created = os.path.exists(DIR_TRIMMED_CLUSTER)
print("Output trimming path created:", out_created, file=sys.stderr)

################################################################################
#Setup out files
################################################################################

if num_chunks == 1:
    # Single chunk => chunk wildcard will be an empty string
    CHUNK_LIST = [""]
else:
    # Multi-chunk => chunk wildcard has suffixes .part_001, .part_002, etc.
    CHUNK_LIST = [f".part_{i:03}" for i in range(0, num_chunks)]


#get all samples from fastq Directory using the fastq2json.py scripts, then just
#load the json file with the samples
FILES = json.load(open(samples))
READS = ["R1", "R2"]  # Paired-end reads
ALL_FASTQ = []
pattern = re.compile(r'[ ._]')

# Normalize sample names by removing underscores
normalized_keys = {}
NEW_FILES = {}
for key in FILES.keys():
    new_key = key.replace("_", "")  # Remove underscores
    if new_key in normalized_keys:
        print(f"Error: Sample name conflict after removing underscores: '{key}' and '{normalized_keys[new_key]}' both map to '{new_key}'", file=sys.stderr)
        sys.exit(1)
    normalized_keys[new_key] = key
    NEW_FILES[new_key] = FILES[key]
FILES = NEW_FILES
ALL_SAMPLES = sorted(FILES.keys())

if fastqs_dir:
    for key in FILES.keys():
        FILES[key]["R1"] = [os.path.abspath(os.path.join(fastqs_dir, i)) for i in FILES[key]["R1"]]
        FILES[key]["R2"] = [os.path.abspath(os.path.join(fastqs_dir, i)) for i in FILES[key]["R2"]]
        
invalid_keys = [key for key in FILES.keys() if pattern.search(key)]
if invalid_keys:
    print("Error: The following JSON file sample names contain invalid characters (spaces, periods, or underscores):\n" + "\n".join(invalid_keys), file=sys.stderr)
    sys.exit()

# Subset FASTQ files if requested by the user (note: these files aren't removed after processing)
if temp_dir and subset_n:
    os.makedirs(temp_dir, exist_ok=True)
    SUBSET_FILES = {}

    def subset_fastq(input_file, output_file, n, buffer_size=100000):
        """Subset a FASTQ file to the first n reads and write to output efficiently."""
        is_gz = input_file.endswith('.gz')
        open_in = gzip.open if is_gz else open
        open_out = lambda f, mode: gzip.open(f, mode, compresslevel=1) if f.endswith('.gz') else open(f, mode)
        with open_in(input_file, 'rt') as infile, open_out(output_file, 'wt') as outfile:
            buffer = []
            read_count = 0
            while read_count < n:
                lines = [infile.readline() for _ in range(4)]  # Read one full read (4 lines)
                if not lines[0]:  # End of file
                    break
                buffer.extend(lines)
                read_count += 1
                if len(buffer) >= buffer_size:  # Write in chunks to optimize I/O
                    outfile.writelines(buffer)
                    buffer = []
            if buffer:  # Write any remaining data in the buffer
                outfile.writelines(buffer)

    for SAMPLE, file in FILES.items():
        SUBSET_FILES[SAMPLE] = {"R1": [], "R2": []}
        for read_type in READS:
            for fastq_file in file.get(read_type, []):
                subset_fastq_path = os.path.join(temp_dir, "subset" + os.path.basename(fastq_file))
                subset_fastq(fastq_file, subset_fastq_path, subset_n)
                SUBSET_FILES[SAMPLE][read_type].append(subset_fastq_path)

    FILES = SUBSET_FILES  # Replace FILES with the subsetted versions

        
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend(
        [i for i in file.get('R1')]
    )
    ALL_FASTQ.extend(
        [i for i in file.get('R2')]
    )


def get_fileids(sample):
    """Numbered indices for R1 files in that sample (should be the same as R2)."""
    return range(len(FILES[sample]["R1"]))

def get_input_r1(wc):
    """
    If wc.chunk == "", return the original R1 file from FILES.
    If wc.chunk != "", return a chunk in split_fastq/...
    """
    # e.g. original path = data/sampleA_R1.fastq.gz
    original = FILES[wc.sample]["R1"][int(wc.fileid)]

    if wc.chunk == "":
        return original
    else:
        # chunked file path
        return f"split_fastq/{wc.sample}_{wc.fileid}_R1{wc.chunk}.fastq.gz"

def get_input_r2(wc):
    """
    Same logic for R2.
    """
    original = FILES[wc.sample]["R2"][int(wc.fileid)]

    if wc.chunk == "":
        return original
    else:
        return f"split_fastq/{wc.sample}_{wc.fileid}_R2{wc.chunk}.fastq.gz"


# Check if any string in the list of samples has a space
if any(" " in string for string in ALL_SAMPLES):
    print("Sample names in json file cannot contain spaces", file=sys.stderr)
    sys.exit()

def get_file_count(sample):
    return len(FILES[sample]["R1"])

def chunked_fastq_path(sample, fileid, read, chunk):
    """
    Final path of chunked FASTQ: e.g. workup/split_fastq/sample_0_R1.part_001.fastq.gz
    If chunk="" => workup/split_fastq/sample_0_R1.fastq.gz (no .part_001 suffix).
    """
    base = f"{sample}_{fileid}_{read}{chunk}.fastq.gz"
    return os.path.join(DIR_WORKUP, "split_fastq", base)

# We create expansions for the final pigz-compressed chunked fastqs for R1 + R2:
SPLIT_FASTQ_GZ = []
for s in ALL_SAMPLES:
    for i in range(get_file_count(s)):
        # R1 expansions
        for c in CHUNK_LIST:
            SPLIT_FASTQ_GZ.append(chunked_fastq_path(s, i, "R1", c))
    for i in range(get_file_count(s)):
        # R2 expansions
        for c in CHUNK_LIST:
            SPLIT_FASTQ_GZ.append(chunked_fastq_path(s, i, "R2", c))
if num_chunks == 1 or not keep_fastq_chunks:
    SPLIT_FASTQ_GZ = []


CONFIG = [os.path.join(DIR_LOGS, "config_" + run_date + "yaml")]

TRIM = [
    os.path.join(DIR_WORKUP, "trimmed", f"{sample}_{fileid}_{read}{chunk}.trimmed.fastq.gz")
    for sample in ALL_SAMPLES
    for fileid in get_fileids(sample)
    for read in READS
    for chunk in CHUNK_LIST
]

# Function to standardize the output file base name
def standardize_name(fq):
    base = os.path.basename(fq)
    if base.endswith('.gz'):
        base = os.path.splitext(base)[0]  # Strip .gz
    base = os.path.splitext(base)[0]  # Strip .fastq, .fq, etc.
    return base

outputs_FASTQC = {fq: standardize_name(fq) for fq in ALL_FASTQ}
FASTQC = expand(os.path.join(DIR_WORKUP, "qc", "{sample}_fastqc.html"), sample=outputs_FASTQC.values())


FASTQC_POST_TRIM = [
    os.path.join(DIR_WORKUP, "qc", "post_trim", f"{sample}_{fileid}_R1{chunk}.trimmed_fastqc.html")
    for sample in ALL_SAMPLES
    for fileid in get_fileids(sample)
    for read in READS
    for chunk in CHUNK_LIST
]



BARCODEID = []
BARCODEID.extend([
    os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R1{chunk}.fastq.gz")
    for sample in ALL_SAMPLES
    for fileid in get_fileids(sample)
    for chunk in CHUNK_LIST
])
BARCODEID.extend([
    os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R2{chunk}.fastq.gz")
    for sample in ALL_SAMPLES
    for fileid in get_fileids(sample)
    for chunk in CHUNK_LIST
])


BARCODEID_MAP = [
    os.path.join(DIR_WORKUP, "assigned", f"{sample}_mapping.txt")                
    for sample in ALL_SAMPLES                               
]

LIGEFF = [
    os.path.join(DIR_WORKUP, f"{sample}.ligation_efficiency.txt")
    for sample in ALL_SAMPLES
]

if not tag_prefix_order:
    LIGEFF = []


EXCLUSION_LOG = []
if use_exclusion_filter:
    EXCLUSION_LOG.extend([
        os.path.join(DIR_WORKUP, "exclusion", f"{sample}_{fileid}.exclusion_alignments{chunk}.log")
        for sample in ALL_SAMPLES
        for fileid in get_fileids(sample)
        for chunk in CHUNK_LIST
    ])
    EXCLUSION_LOG.extend([
        os.path.join(DIR_WORKUP, f"{sample}.exclusion_stats.txt")
        for sample in ALL_SAMPLES
    ])


OUT_FASTQS = [
    os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R1{chunk}.fastq.gz")
    for sample in ALL_SAMPLES                      
    for fileid in get_fileids(sample)
    for chunk in CHUNK_LIST                                    
]
OUT_FASTQS.extend([              
    os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R2{chunk}.fastq.gz")
    for sample in ALL_SAMPLES
    for fileid in get_fileids(sample)
    for chunk in CHUNK_LIST
])                                                       

OUT_STAR = expand(
    os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Aligned.sortedByCoord.out.bam"),
    sample=ALL_SAMPLES
)
OUT_STAR.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.tsv"), sample=ALL_SAMPLES))                
OUT_STAR.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.tsv"), sample=ALL_SAMPLES))                
OUT_STAR.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "SJ", "raw", "barcodes.tsv"), sample=ALL_SAMPLES))                

OUT_STAR_BARCODE_MAPPING = []

OUT_STAR_BARCODE_MAPPING.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.ids.txt"), sample=ALL_SAMPLES))
OUT_STAR_BARCODE_MAPPING.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.ids.txt"), sample=ALL_SAMPLES))
OUT_STAR_BARCODE_MAPPING.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "SJ", "raw", "barcodes.ids.txt"), sample=ALL_SAMPLES))


if not use_star:
    OUT_STAR = []
    OUT_STAR_BARCODE_MAPPING = []

if use_star:
    if split_allele_analysis:
        OUT_STAR.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_), sample=ALL_SAMPLES))
        OUT_STAR.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2), sample=ALL_SAMPLES))
    OUT_ANNDATA = expand(
        [os.path.join(
            DIR_WORKUP,
            "results",
            "anndatas",
            "{sample}_anndata",
            "starsolo_" + align_id,
            "Gene",
            adata_filename),
         os.path.join(
            DIR_WORKUP,
            "results",
            "anndatas",
            "{sample}_anndata",
            "starsolo_" + align_id,
            "GeneFull",
            adata_filename)],
        sample=ALL_SAMPLES
    )
    if merge_anndata_files:
        OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "Gene", adata_filename)])
        OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "GeneFull", adata_filename)])
        OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "Gene", "sample_order.txt")])
        OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "GeneFull", "sample_order.txt")])
    if split_allele_analysis:
        OUT_ANNDATA.extend(expand(
            [os.path.join(
                DIR_WORKUP,
                "results",
                "anndatas",
                "{sample}_anndata",
                "starsolo_" + align_id,
                vcf_s1_,
                "Gene",
                adata_filename),
             os.path.join(
                DIR_WORKUP,
                "results",
                "anndatas",
                "{sample}_anndata",
                "starsolo_" + align_id,
                vcf_s1_,
                "GeneFull",
                adata_filename)],
            sample=ALL_SAMPLES
        ))
        OUT_ANNDATA.extend(expand(
            [os.path.join(
                DIR_WORKUP,
                "results",
                "anndatas",
                "{sample}_anndata",
                "starsolo_" + align_id,
                vcf_s2,
                "Gene",
                adata_filename),
             os.path.join(
                DIR_WORKUP,
                "results",
                "anndatas",
                "{sample}_anndata",
                "starsolo_" + align_id,
                vcf_s2,
                "GeneFull",
                adata_filename)],
            sample=ALL_SAMPLES
        ))
        if merge_anndata_files:
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "Gene", adata_filename)])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "GeneFull", adata_filename)])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "Gene", adata_filename)])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "GeneFull", adata_filename)])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "Gene", "sample_order.txt")])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "GeneFull", "sample_order.txt")])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "Gene", "sample_order.txt")])
            OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "GeneFull", "sample_order.txt")])

else:
    OUT_ANNDATA = []


STAR_GENOME_FASTA = []
STAR_GENOME_GTF = []
if use_star and genome_fasta and genome_gtf and not star_index:
    STAR_GENOME_FASTA = [genome_fasta]
    STAR_GENOME_GTF = [genome_gtf]

# features files:

OUT_FEATURES = []
if use_features:
    OUT_FEATURES = expand(
        os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.mtx"),
        sample=ALL_SAMPLES
    )
    OUT_FEATURES.extend(expand(os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.barcodes.ids.txt"), sample=ALL_SAMPLES))
    OUT_ANNDATA.extend(expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "features", adata_filename), sample=ALL_SAMPLES))
    if merge_anndata_files:
        OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "features", adata_filename)])
        OUT_ANNDATA.extend([os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "features", "sample_order.txt")])


# kallisto files:

OUT_KALLISTO = expand(
    os.path.join(DIR_WORKUP, "results", "output_{sample}_kallisto_" + align_id, "counts_unfiltered"),
    sample=ALL_SAMPLES
)

if not use_kallisto:
    OUT_KALLISTO = []

################################################################################
#Functions for formatting of input file names and output file names
################################################################################


#returns path to supplied pre-generated index or a newly created index
star_index_dir_new = os.path.join(DIR_WORKUP, "index", "star", "genome")
def get_star_index(wildcards):
    if star_index:
        return star_index
    else:
        return star_index_dir_new


#kallisto
kallisto_index_dir_new = os.path.join(DIR_WORKUP, "index", "kallisto")
kallisto_index_new = os.path.join(kallisto_index_dir_new, "index.idx")
kallisto_t2g_new = os.path.join(kallisto_index_dir_new, "t2g.txt")
kallisto_c1_new = os.path.join(kallisto_index_dir_new, "cdna.txt")
kallisto_c2_new = os.path.join(kallisto_index_dir_new, "nascent.txt")
kallisto_f1_new = os.path.join(kallisto_index_dir_new, "cdna.fa")
kallisto_f2_new = os.path.join(kallisto_index_dir_new, "nascent.fa")
def get_kallisto_index(wildcards):
    if kallisto_index:
        return kallisto_index
    else:
        return kallisto_index_new
def get_kallisto_t2g(wildcards):
    if kallisto_index:
        return kallisto_index_t2g
    else:
        return kallisto_t2g_new
def get_kallisto_c1(wildcards):
    if kallisto_index:
        return kallisto_index_c1
    else:
        return kallisto_c1_new
def get_kallisto_c2(wildcards):
    if kallisto_index:
        return kallisto_index_c2
    else:
        return kallisto_c2_new


#returns path to VCF file
def get_vcf_file(wildcards):
    if not vcf_file:
        return []
    if use_vcf_wasp_workflow and vcf_s1 and vcf_s2:
        return [os.path.join(DIR_WORKUP, "hybrid.vcf")]
    return [vcf_file]

#variable describing how the VCF file should be inputted
vcf_file_input_param = "\"" + os.path.join(DIR_WORKUP, "hybrid.vcf") + "\"" if vcf_s1 and vcf_s2 and use_vcf_wasp_workflow else "<(zcat -f \"" + vcf_file + "\")"


#retrieves the fully preprocessed paired-end reads files to be used in read mapping
def get_mapping_input_r1(wildcards):
    """
    Returns a list of all R1 FASTQ files for the given sample across all fileids and chunks.
    """
    sample = wildcards.sample
    r1_files = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            if not use_exclusion_filter:
                filepath = os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R1{chunk}.fastq.gz")
            else:
                filepath = os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R1{chunk}.filtered.fastq.gz")
            r1_files.append(filepath)
    return r1_files

def get_mapping_input_r2(wildcards):
    """
    Returns a list of all R2 FASTQ files for the given sample across all fileids and chunks.
    """
    sample = wildcards.sample
    r2_files = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            if not use_exclusion_filter:
                filepath = os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R2{chunk}.fastq.gz")
            else:
                filepath = os.path.join(DIR_WORKUP, "fastqs", f"{sample}_{fileid}_R2{chunk}.filtered.fastq.gz")
            r2_files.append(filepath)
    return r2_files

def get_mapping_input_names(wildcards):
    """
    Returns a list of all "chunk names" for a given pair of read files.
    """
    sample = wildcards.sample
    chunk_names = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            chunk_name = f"ID:{sample}" # f"ID:{sample}\tXF:{fileid}\tXP:{chunk.lstrip('.')}"
            chunk_names.append(chunk_name)
    return chunk_names
  
def get_group_from_input_r2_name(input_string):
    """
    Returns a list of all "chunk names" for a given R2 read file.
    """
    match = re.match(r"(?P<sample>[^_]+)_(?P<fileid>[^_]+)_R2(?P<chunk>\.[^.]+)?(\..+)?", os.path.basename(input_string))
    if match:
        sample = match.group("sample")
        fileid = match.group("fileid")
        chunk = match.group("chunk")
        return f"ID:{sample}"  # f"ID:{sample}\tXF:{fileid}\tXP:{chunk.lstrip('.')}"
    else:
        return "ID:undefined"


#retrieves the fully preprocessed paired-end reads files to be used for split allele analysis
def get_mapping_input_r1_split1(wildcards):
    """
    Returns a list of all R1 FASTQ files for the given sample across all fileids and chunks.
    """
    sample = wildcards.sample
    r1_files = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            filepath = os.path.join(DIR_WORKUP, "chunks", f"{sample}_{fileid}_R1{chunk}.split1.fastq.gz")
            r1_files.append(filepath)
    return r1_files

def get_mapping_input_r2_split1(wildcards):
    """
    Returns a list of all R2 FASTQ files for the given sample across all fileids and chunks.
    """
    sample = wildcards.sample
    r2_files = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            filepath = os.path.join(DIR_WORKUP, "chunks", f"{sample}_{fileid}_R2{chunk}.split1.fastq.gz")
            r2_files.append(filepath)
    return r2_files
  

def get_mapping_input_r1_split2(wildcards):
    """
    Returns a list of all R1 FASTQ files for the given sample across all fileids and chunks.
    """
    sample = wildcards.sample
    r1_files = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            filepath = os.path.join(DIR_WORKUP, "chunks", f"{sample}_{fileid}_R1{chunk}.split2.fastq.gz")
            r1_files.append(filepath)
    return r1_files

def get_mapping_input_r2_split2(wildcards):
    """
    Returns a list of all R2 FASTQ files for the given sample across all fileids and chunks.
    """
    sample = wildcards.sample
    r2_files = []
    for fileid in get_fileids(sample):
        for chunk in CHUNK_LIST:
            filepath = os.path.join(DIR_WORKUP, "chunks", f"{sample}_{fileid}_R2{chunk}.split2.fastq.gz")
            r2_files.append(filepath)
    return r2_files

#formats path to only consist of the filename
def get_basename():
    return splitext(basename(file_path))[0]

################################################################################
################################################################################
#Rule all
################################################################################
################################################################################


rule all:
    input: CONFIG + FASTQC + FASTQC_POST_TRIM + BARCODEID_MAP + SPLIT_FASTQ_GZ + EXCLUSION_LOG + \
           LIGEFF + OUT_FASTQS + OUT_STAR + OUT_ANNDATA + OUT_KALLISTO + OUT_FEATURES + OUT_STAR_BARCODE_MAPPING

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')


################################################################################
#Trimming and barcode identification
################################################################################



rule log_config:
    '''Copy config.yaml and place in logs folder with the date run
    '''
    input:
        config_path
    output:
        os.path.join(DIR_LOGS, "config_" + run_date + "yaml")
    shell:
        '''
        cp "{input}" "{output}"
        '''


if num_chunks > 1:
    rule split_fastq_into_parts_r1:
        """
        Splits original FASTQ (R1) into chunks.
        """
        input:
            r1=lambda wc: FILES[wc.sample]['R1'][int(wc.fileid)]
        output:
            temp(expand([os.path.join(DIR_WORKUP, "split_fastq", "{{sample}}_{{fileid}}_R1{splitid}.fastq")], splitid=CHUNK_LIST))
        threads:
            4
        priority:
            8
        params:
            dir=os.path.join(DIR_WORKUP, "split_fastq")
        log:
            os.path.join(DIR_LOGS, "{sample}_{fileid}_split_r1_fastq.log")
        shell:
            """
            bash "{split_fastq}" "{input.r1}" {num_chunks} "{params.dir}" "{wildcards.sample}_{wildcards.fileid}_R1.part_" {threads} &> "{log}"
            """
            

    rule split_fastq_into_parts_r2:
        """
        Splits original FASTQ (R2) into chunks.
        """
        input:
            r2=lambda wc: FILES[wc.sample]['R2'][int(wc.fileid)]
        output:
            temp(expand([os.path.join(DIR_WORKUP, "split_fastq", "{{sample}}_{{fileid}}_R2{splitid}.fastq")], splitid=CHUNK_LIST))
        threads:
            4
        priority:
            8
        params:
            dir=os.path.join(DIR_WORKUP, "split_fastq")
        log:
            os.path.join(DIR_LOGS, "{sample}_{fileid}_split_r2_fastq.log")
        shell:
            """
            bash "{split_fastq}" "{input.r2}" {num_chunks} "{params.dir}" "{wildcards.sample}_{wildcards.fileid}_R2.part_" {threads} &> "{log}"
            """

    rule compress_split_fastq:
        input:
            r1 = os.path.join(DIR_WORKUP, "split_fastq", "{sample}_{fileid}_R1{chunk}.fastq"),
            r2 = os.path.join(DIR_WORKUP, "split_fastq", "{sample}_{fileid}_R2{chunk}.fastq")
        output:
            r1 = os.path.join(DIR_WORKUP, "split_fastq", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
            r2 = os.path.join(DIR_WORKUP, "split_fastq", "{sample}_{fileid}_R2{chunk}.fastq.gz")
        conda:
            conda_env
        threads:
            4*thread_factor
        log:
            os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}pigz.log")
        shell:
            '''
            pigz -k -p {threads} "{input.r1}" &> "{log}"
            pigz -k -p {threads} "{input.r2}" &>> "{log}"
            '''


rule initial_fastqc:
    '''Initial QC of raw input reads
    '''
    input:
        lambda wildcards: [fq for fq, base in outputs_FASTQC.items() if base == wildcards.sample][0]
    output:
        html=os.path.join(DIR_WORKUP, "qc", "{sample}_fastqc.html")
    threads:
        1
    log:
        os.path.join(DIR_LOGS, "{sample}.initial_fastqc.log")
    params:
        dir_qc = os.path.join(DIR_WORKUP, "qc")
    shell:
        '''
        fastqc "{input}" -o "{params.dir_qc}" &> "{log}"
        '''



rule adaptor_trimming_pe:
   '''Initial trimming of index adapters from raw input reads for chunked data'''
    input:
        reads1 = os.path.join(DIR_WORKUP, "split_fastq", "{sample}_{fileid}_R1{chunk}.fastq") if num_chunks > 1 else lambda wc: FILES[wc.sample]['R1'][int(wc.fileid)],
        reads2 = os.path.join(DIR_WORKUP, "split_fastq", "{sample}_{fileid}_R2{chunk}.fastq") if num_chunks > 1 else lambda wc: FILES[wc.sample]['R2'][int(wc.fileid)]
    output:
        r1 = temp(os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}_R1{chunk}.trimmed.fastq.gz")),
        r2 = temp(os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}_R2{chunk}.trimmed.fastq.gz")),
        qc = os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}{chunk}.trimmed.qc.txt")
   threads:
        4*thread_factor
   priority:
        7
   log:
        os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.cutadapt_adaptor_trimming.log")
   conda:
        conda_env
   shell:
       '''
        (cutadapt \
       -j {threads} \
       -m 31 \
       --nextseq-trim=20 \
       -a {R1_adapter} \
       -A {R2_adapter} \
       -o "{output.r1}" -p "{output.r2}" \
        "{input.reads1}" "{input.reads2}" > "{output.qc}") &> "{log}"
        '''



rule post_trim_fastqc:
    '''Post-adapter-trimming QC of trimmed reads
    '''
    input:
        r1 = os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}_R1{chunk}.trimmed.fastq.gz"),
        r2 = os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}_R2{chunk}.trimmed.fastq.gz"),
    output:
        html1=os.path.join(DIR_WORKUP, "qc", "post_trim", "{sample}_{fileid}_R1{chunk}.trimmed_fastqc.html"),
        html2=os.path.join(DIR_WORKUP, "qc", "post_trim", "{sample}_{fileid}_R2{chunk}.trimmed_fastqc.html")
    threads:
        2
    log:
        os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.post_trim_fastqc.log")
    params:
        dir_qc = os.path.join(DIR_WORKUP, "qc", "post_trim")
    conda:
        conda_env
    shell:
        '''
        fastqc -t {threads} -o "{params.dir_qc}" "{input.r1}" "{input.r2}" &> "{log}"
        '''


# Set up the UMI extraction pattern
umi_extraction_pattern = umi_extraction_pattern.replace("{","{{").replace("}","}}").replace("R1","{params.r1_name}").replace("R2","{params.r2_name}")


rule splitcode_barcodeID:
    '''splitcode for barcode identification/assignment/formatting
       Output:
       R1: 20-bp barcode + 10-bp UMI + sequence-to-be-mapped
       R2: sequence-to-be-mapped
    '''
    input:
        r1 = os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}_R1{chunk}.trimmed.fastq.gz"),
        r2 = os.path.join(DIR_WORKUP, "trimmed", "{sample}_{fileid}_R2{chunk}.trimmed.fastq.gz")
    output:
        r1 = os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
        r2 = os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.fastq.gz"),
        u1 = temp(os.path.join(DIR_WORKUP, "assigned", "{sample}_{fileid}_unassigned_R1{chunk}.fq.gz")),
        map = temp(os.path.join(DIR_WORKUP, "assigned", "{sample}_{fileid}_mapping{chunk}.txt"))
    threads:
        4*thread_factor
    log:
        os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.splitcode_barcodeID.log")
    params:
        r1_name=lambda wc: os.path.join(DIR_WORKUP, "fastqs", f"{wc.sample}_{wc.fileid}_R1{wc.chunk}"),
        r2_name=lambda wc: os.path.join(DIR_WORKUP, "fastqs", f"{wc.sample}_{wc.fileid}_R2{wc.chunk}")
    conda:
        conda_env
    shell:
        '''
        splitcode \
        --assign \
        --x-only \
        --nFastqs=2 \
        --empty N \
        -x "{}" \
        --gzip \
        --mod-names \
        --bclen=20 \
        -t {{threads}} \
        --unassigned="{{output.u1}}," \
        -c "{{splitcode_bID}}" \
        --mapping="{{output.map}}" \
        "{{input.r1}}" "{{input.r2}}" &> "{{log}}"
        '''.format(umi_extraction_pattern)



rule calc_ligation_efficiency:
    '''Calculate ligation efficiency
       Note: We get the expected barcode tag order by looking
       at the line after @keep-grp in the splitcode config file
    '''
    input:
        barcodes = os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
        u1 = os.path.join(DIR_WORKUP, "assigned", "{sample}_{fileid}_unassigned_R1{chunk}.fq.gz")
    output:
        temp(os.path.join(DIR_WORKUP, "chunks", "{sample}_{fileid}.ligation_efficiency{chunk}.txt"))
    log:
        os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.calc_ligation_efficiency.log")
    threads:
        1
    shell:
        '''
        (bash "{ligeff}" \
        "{input.barcodes}" "{input.u1}" \
        "{tag_prefix_order}" \
        > "{output}") &> "{log}"
        '''


rule concatenate_ligeff_results:
    """
    Concatenates all ligation efficiency files across fileids and chunks for each sample.
    """
    input:
        # For each sample, gather all {fileid} and {chunk} combinations
        ligeff = lambda wc: [
            os.path.join(DIR_WORKUP, "chunks", f"{wc.sample}_{fileid}.ligation_efficiency{chunk}.txt")
            for fileid in get_fileids(wc.sample)
            for chunk in CHUNK_LIST
        ]
    output:
        ligeff=os.path.join(DIR_WORKUP, "{sample}.ligation_efficiency.txt")
    log:
        os.path.join(DIR_LOGS, "{sample}.concat_ligeff.log")
    shell:
        '''
        python "{merge_ligeff}" -o "{output.ligeff}" {input.ligeff} &> "{log}"
        '''

rule concatenate_exclusion_results:
    """
    Concatenates all exclusion stats files across fileids and chunks for each sample.
    """
    input:
        # For each sample, gather all {fileid} and {chunk} combinations
        stats = lambda wc: sum([
            [os.path.join(DIR_WORKUP, "exclusion", f"{wc.sample}_{fileid}.exclusion_alignments{chunk}.log"), os.path.join(DIR_WORKUP, "exclusion", f"{wc.sample}_{fileid}.exclusion_numreads{chunk}.log")]
            for fileid in get_fileids(wc.sample)
            for chunk in CHUNK_LIST
        ], [])
    output:
        stats=os.path.join(DIR_WORKUP, "{sample}.exclusion_stats.txt")
    log:
        os.path.join(DIR_LOGS, "{sample}.concat_exclusion_stats.log")
    shell:
        '''
        python "{merge_exclusion_stats}" "{output.stats}" {input.stats} &> "{log}"
        '''


rule set_barcode_mapping_file:
    """
    Sets the final barcode mapping file (since multiple might exist if there are multiple chunks even though they're identical).
    """
    input:
        # For each sample, gather all {fileid} and {chunk} combinations
        map = lambda wc: [                                         
            os.path.join(DIR_WORKUP, "assigned", f"{wc.sample}_{fileid}_mapping{chunk}.txt")
            for fileid in get_fileids(wc.sample)                                               
            for chunk in CHUNK_LIST                                              
        ]
    output:
        map=os.path.join(DIR_WORKUP, "assigned", "{sample}_mapping.txt")
    log:
        os.path.join(DIR_LOGS, "{sample}.set_mapping_file.log")
    shell:
        '''
        (bash -c 'cat "$1"' _ {input.map} > "{output.map}") &> "{log}"
        '''


if num_chunks == 1 or not chunk_star_alignment:
    
    rule starsolo_align:
        '''STARsolo to paired-end-align FASTQ files and generate count matrices and BAM files
           We ignore the first 30-bp (20 bp barcode + 10 bp UMI) of R1 via clip5pNbases
        '''
        input:
            r1=get_mapping_input_r1,
            r2=get_mapping_input_r2,
            index=get_star_index,
            vcf=get_vcf_file
        output:
            bam=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Aligned.sortedByCoord.out.bam"),
            gene=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.tsv"),
            genefull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.tsv"),
            sj=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "SJ", "raw", "barcodes.tsv")
        threads:
            18
        log:
            os.path.join(DIR_LOGS, "{sample}.starsolo." + align_id + ".log")
        params:
            wasp="--waspOutputMode SAMtag --varVCFfile " + vcf_file_input_param if use_vcf_wasp_workflow else "",
            wasp_sam_attr="NH HI AS nM NM NM MD vA vG vW" if use_vcf_wasp_workflow else "",
            id="{sample}",
            additional=star_additional_params,
            out=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id),
            r1=lambda wildcards, input: ",".join([f'"{x}"' for x in input.r1]),
            r2=lambda wildcards, input: ",".join([f'"{x}"' for x in input.r2]),
            nthreads=lambda wildcards, input, threads: threads - 1  # Since samtools uses an additional thread as the main thread
        conda:
            conda_env
        shell:
            '''
            STAR {params.wasp} {params.additional} \
            --genomeDir "{input.index}" \
            --outSAMattributes CR UR {params.wasp_sam_attr} \
            --soloCellReadStats Standard \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist None \
            --soloBarcodeMate 1 \
            --clip5pNbases 30 0 \
            --soloCBstart 1 --soloCBlen 20 \
            --soloUMIstart 21 --soloUMIlen 10 \
            --soloUMIdedup Exact --soloUMIfiltering MultiGeneUMI_All \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --soloFeatures Gene GeneFull SJ \
            --outSAMattrRGline "ID:{params.id}" "SM:{params.id}" \
            --outFileNamePrefix "{params.out}/" \
            --readFilesIn {params.r1} {params.r2} &> "{log}"
            samtools index -@ {params.nthreads} "{output.bam}" &>> "{log}"
            '''

else:
  
    rule starsolo_generate_manifest:
        input:
            r1=get_mapping_input_r1,
            r2=get_mapping_input_r2
        output:
            os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "manifest.tsv")
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "{sample}.manifest.starsolo." + align_id + ".log")
        params:
            id="{sample}",
            r1=lambda wildcards, input: ",".join(input.r1),
            r2=lambda wildcards, input: ",".join(input.r2),
            rg=lambda wildcards, input: ",".join([get_group_from_input_r2_name(s) for s in input.r2])
        conda:
            conda_env
        shell:
            '''
            (paste \
            <(echo "{params.r1}" | tr ',' '\n') \
            <(echo "{params.r2}" | tr ',' '\n') \
            <(echo "{params.rg}" | tr ',' '\n') > "{output}") &> "{log}"
            '''
            
    rule starsolo_align_manifest:
        '''STARsolo to paired-end-align FASTQ files and generate count matrices and BAM files
           We ignore the first 30-bp (20 bp barcode + 10 bp UMI) of R1 via clip5pNbases
           Individual chunks are aligned; not the entire dataset.
        '''
        input:
            r1=get_mapping_input_r1,
            r2=get_mapping_input_r2,
            manifest=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "manifest.tsv"),
            index=get_star_index,
            vcf=get_vcf_file
        output:
            bam=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Aligned.sortedByCoord.out.bam"),
            gene=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.tsv"),
            genefull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.tsv"),
            sj=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "SJ", "raw", "barcodes.tsv")
        threads:
            18*thread_factor                                             
        log:    
            os.path.join(DIR_LOGS, "{sample}.starsolo." + align_id + ".log") 
        params:   
            wasp="--waspOutputMode SAMtag --varVCFfile " + vcf_file_input_param if use_vcf_wasp_workflow else "",
            wasp_sam_attr="NH HI AS nM NM NM MD vA vG vW" if use_vcf_wasp_workflow else "",
            id="{sample}",
            additional=star_additional_params,
            out=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id),
            nthreads=lambda wildcards, input, threads: threads - 1  # Since samtools uses an additional thread as the main thread
        conda: 
            conda_env
        shell:     
            '''
            STAR {params.wasp} {params.additional} \
            --genomeDir "{input.index}" \
            --outSAMattributes CR UR {params.wasp_sam_attr} \
            --soloCellReadStats Standard \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist None \
            --soloBarcodeMate 1 \
            --clip5pNbases 30 0 \
            --soloCBstart 1 --soloCBlen 20 \
            --soloUMIstart 21 --soloUMIlen 10 \
            --soloUMIdedup Exact --soloUMIfiltering MultiGeneUMI_All \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --soloFeatures Gene GeneFull SJ \
            --outSAMattrRGline "ID:{params.id}" "SM:{params.id}" \
            --outFileNamePrefix "{params.out}/" \
            --readFilesManifest "{input.manifest}" &> "{log}"
            samtools index -@ {params.nthreads} "{output.bam}" &>> "{log}"
            '''

rule decode_barcodes_starsolo:
    '''Map/decode barcodes from STARsolo to their barcode IDs
    '''
    input:
        map=os.path.join(DIR_WORKUP, "assigned", "{sample}_mapping.txt"),
        gene=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.tsv"),
        genefull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.tsv"),
        sj=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "SJ", "raw", "barcodes.tsv")
    output:
        gene=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.ids.txt"),
        genefull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.ids.txt"),
        sj=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "SJ", "raw", "barcodes.ids.txt")
    threads:
        1
    log:
        os.path.join(DIR_LOGS, "{sample}.decode_barcode.starsolo." + align_id + ".log")
    conda:
        conda_env
    shell:
        '''
        python "{decode_barcode}" "{input.map}" "{input.gene}" "{output.gene}" &> "{log}"
        python "{decode_barcode}" "{input.map}" "{input.genefull}" "{output.genefull}" &>> "{log}"
        python "{decode_barcode}" "{input.map}" "{input.sj}" "{output.sj}" &>> "{log}"
        '''


rule generate_anndata_starsolo:
    '''Create anndata object from STARsolo matrices
    '''
    input:
        bam=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Aligned.sortedByCoord.out.bam"),
        mapgene=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.ids.txt"),
        mapgenefull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.ids.txt")
    output:
        gene=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, "Gene", adata_filename),
        genefull=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, "GeneFull", adata_filename)
    threads:
        1
    log:
        os.path.join(DIR_LOGS, "{sample}.anndata.starsolo." + align_id + ".log")
    params:
        dir=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw"),
        out=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, "Gene"),
        dirfull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw"),
        outfull=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, "GeneFull")
    conda:
        conda_env
    shell:
        '''
        python "{prep_anndata}" STAR "{params.dir}" "{input.mapgene}" "{params.out}" "" "{select_tags}" &> "{log}"
        python "{prep_anndata}" STAR "{params.dirfull}" "{input.mapgenefull}" "{params.outfull}" "" "{select_tags}" &>> "{log}"
        '''


if merge_anndata_files:
    rule merge_anndata:
        '''Create a merged anndata file
        '''
        input:
            gene=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, "Gene", adata_filename), sample=ALL_SAMPLES),
            genefull=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, "GeneFull", adata_filename), sample=ALL_SAMPLES)
        output:
            gene=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "Gene", adata_filename),
            genefull=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "GeneFull", adata_filename),
            info_gene=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "Gene", "sample_order.txt"),
            info_genefull=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, "GeneFull", "sample_order.txt")
        log:
            os.path.join(DIR_LOGS, "anndata.merged." + align_id + ".log")
        threads:
            1
        conda:
            conda_env
        shell:
            '''
            python "{merge_anndata}" "{output.gene}" {input.gene} &> "{log}"
            python "{merge_anndata}" "{output.genefull}" {input.genefull} &>> "{log}"
            (echo "{input.gene}" > "{output.info_gene}") &>> "{log}"
            (echo "{input.genefull}" > "{output.info_genefull}") &>> "{log}"
            '''


if use_exclusion_filter:

    rule create_exclusion_index:
        '''Create an exclusion index via bowtie2 based on supplied exclusion FASTA file
        '''
        input:
            exclusion_fasta_file
        output:
            os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.1.bt2"), 
            os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.2.bt2"),
            os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.3.bt2"),
            os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.4.bt2"),
            os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.rev.1.bt2"),
            os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.rev.2.bt2")
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "bowtie2.exclusion_index" + ".log")
        conda:
            conda_env
        params:
            prefix=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index")
        shell:
            '''
            bowtie2-build "{input}" "{params.prefix}" &> "{log}"
            '''

    rule get_exclusion_reads:
        '''bowtie2 alignment to R1+R2 FASTQ files to get excluded read names
           We ignore the first 30-bp (20 bp barcode + 10 bp UMI) of R1
        '''
        input:
            index1=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.1.bt2"),
            index2=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.2.bt2"),
            index3=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.3.bt2"),
            index4=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.4.bt2"),
            index5=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.rev.1.bt2"),
            index6=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index.rev.2.bt2"),
            r1=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
            r2=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.fastq.gz")
        output:
            list=temp(os.path.join(DIR_WORKUP, "bowtie2", "{sample}_{fileid}.exclusion_list{chunk}.txt")),
            stats=os.path.join(DIR_WORKUP, "exclusion", "{sample}_{fileid}.exclusion_alignments{chunk}.log")
        threads:
            4*thread_factor
        priority:
            5
        params:
            prefix=os.path.join(DIR_WORKUP, "bowtie2", "exclusion_index")
        conda:
            conda_env
        shell:
            '''
            set +o pipefail
            (bowtie2 -q -p {threads} \
            --no-unal \
            --local \
            -x "{params.prefix}" \
            -1 <(zcat "{input.r1}"|awk '{{if (NR%4==2 || NR%4==0) $0=substr($0,31); print}}') \
            -2 "{input.r2}" \
            -5 0 | samtools view -S | cut -f1 | uniq > "{output.list}") &> "{output.stats}"
            '''

    rule exclude_reads_r1:
        '''Exclude reads based off a list of read names
        '''
        input:
            list=os.path.join(DIR_WORKUP, "bowtie2", "{sample}_{fileid}.exclusion_list{chunk}.txt"),
            r1=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
            r2=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.fastq.gz")
        output:
            r1=(os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.filtered.fastq.gz"))
        threads:
            2
        log:
            os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.exclusion.exclude_reads.log")
        conda:
            conda_env
        shell:
            '''
            (seqkit \
            grep -j {threads} -v -n \
            -f "{input.list}" "{input.r1}" \
            -o "{output.r1}") &> "{log}"
            '''

    rule exclude_reads_r2:
        '''Exclude reads based off a list of read names
        '''
        input:
            list=os.path.join(DIR_WORKUP, "bowtie2", "{sample}_{fileid}.exclusion_list{chunk}.txt"),
            r1=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
            r2=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.fastq.gz")
        output:
            r2=(os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.filtered.fastq.gz")),
            log=temp(os.path.join(DIR_WORKUP, "exclusion", "{sample}_{fileid}.exclusion_numreads{chunk}.log"))
        threads:
            2
        conda:
            conda_env
        shell:
            '''
            (seqkit \
            grep -j {threads} -v -n \
            -f "{input.list}" "{input.r2}" \
            -o "{output.r2}") &> "{output.log}"
            '''

# If both species are non-REF, make a hybrid VCF of the two ALT species provided
if use_vcf_wasp_workflow and vcf_s1 and vcf_s2:

    rule make_hybrid_vcf:
        '''Makes a hybrid non-REF two-species VCF file for use in STAR
        '''
        input:
            vcf=vcf_file
        output:
            vcf=temp(os.path.join(DIR_WORKUP, "hybrid.vcf"))
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "prep_vcf.log")
        conda:
            conda_env
        shell:
            '''
            python "{prep_vcf}" <(zcat -f "{input.vcf}") "{vcf_s1}" "{vcf_s2}" "{output.vcf}" &> "{log}"
            '''

# Split BAM by species
if use_vcf_wasp_workflow and split_allele_analysis:

    rule split_by_allele_1:
        '''Split BAM files into two set of read names based on allele-specific alignment results
        '''
        input:
            bam=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Aligned.sortedByCoord.out.bam")
        output:
            names1=temp(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "split_names_1.txt")),
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "{sample}.starsolo.split1.bam." + align_id + ".log")
        conda:
            conda_env
        shell:
            '''
            set +o pipefail
            (samtools view "{input.bam}" | grep "vW:i:1" | grep -P 'vA:B:c(,1)+(\s|\t|$)' | cut -f1 > "{output.names1}") &> "{log}"
            '''

    rule split_by_allele_2:
        '''Split BAM files into two set of read names based on allele-specific alignment results
        '''
        input:
            bam=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Aligned.sortedByCoord.out.bam")
        output:
            names2=temp(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "split_names_2.txt"))
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "{sample}.starsolo.split2.bam." + align_id + ".log")
        conda:
            conda_env
        shell:
            '''
            set +o pipefail
            (samtools view "{input.bam}" | grep "vW:i:1" | grep -P 'vA:B:c(,2)+(\s|\t|$)' | cut -f1 > "{output.names2}") &> "{log}"
            '''

    rule filter_allelic_reads_1:
        '''Filter FASTQ files for the read names from the allelic splitting
        '''
        input:
            list1=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "split_names_1.txt"),
            r1=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
            r2=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.fastq.gz")
        output:
            r1_split1=temp(os.path.join(DIR_WORKUP, "chunks", "{sample}_{fileid}_R1{chunk}.split1.fastq.gz")),
            r2_split1=temp(os.path.join(DIR_WORKUP, "chunks", "{sample}_{fileid}_R2{chunk}.split1.fastq.gz"))
        threads:
            2
        log:
            os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.filter_allelic_reads_1.log")
        conda:
            conda_env
        shell:
            '''
            (seqkit \
            grep -j {threads} -n \
            -f "{input.list1}" "{input.r1}" \
            -o "{output.r1_split1}" && seqkit \
            grep -j {threads} -n \
            -f "{input.list1}" "{input.r2}" \
            -o "{output.r2_split1}") &> "{log}"
            '''

    rule filter_allelic_reads_2:
        '''Filter FASTQ files for the read names from the allelic splitting
        '''
        input:
            list2=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "split_names_2.txt"),
            r1=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R1{chunk}.fastq.gz"),
            r2=os.path.join(DIR_WORKUP, "fastqs", "{sample}_{fileid}_R2{chunk}.fastq.gz")
        output:
            r1_split2=temp(os.path.join(DIR_WORKUP, "chunks", "{sample}_{fileid}_R1{chunk}.split2.fastq.gz")),
            r2_split2=temp(os.path.join(DIR_WORKUP, "chunks", "{sample}_{fileid}_R2{chunk}.split2.fastq.gz"))
        threads:
            2
        log:
            os.path.join(DIR_LOGS, "{sample}_{fileid}{chunk}.filter_allelic_reads_2.log")
        conda:
            conda_env
        shell:
            '''
            (seqkit \
            grep -j {threads} -n \
            -f "{input.list2}" "{input.r1}" \
            -o "{output.r1_split2}" && seqkit \
            grep -j {threads} -n \
            -f "{input.list2}" "{input.r2}" \
            -o "{output.r2_split2}") &> "{log}"
            '''

    rule starsolo_align_post_splitting_1:
        '''STARsolo to paired-end-align FASTQ files and generate count matrices and BAM files
           We ignore the first 30-bp (20 bp barcode + 10 bp UMI) of R1 via clip5pNbases
           This rule is run on FASTQ files that have been split by allele
        '''
        input:
            r1_split1=get_mapping_input_r1_split1,
            r2_split1=get_mapping_input_r2_split1,
            index=get_star_index,
            vcf=get_vcf_file
        output:
            s1=directory(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_))
        threads:
            16
        log:
            os.path.join(DIR_LOGS, "{sample}.post_split1.starsolo." + align_id + ".log")
        params:
            wasp="--waspOutputMode SAMtag --varVCFfile " + vcf_file_input_param if use_vcf_wasp_workflow else "",
            wasp_sam_attr="NH HI AS nM NM NM MD vA vG vW" if use_vcf_wasp_workflow else "",
            id="{sample}",
            additional=star_additional_params,
            b1=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_, "Aligned.sortedByCoord.out.bam"),
            b2=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2, "Aligned.sortedByCoord.out.bam"),
            r1=lambda wildcards, input: ",".join([f'"{x}"' for x in input.r1_split1]),
            r2=lambda wildcards, input: ",".join([f'"{x}"' for x in input.r2_split1]),
            nthreads=lambda wildcards, input, threads: threads - 1  # Since samtools uses an additional thread as the main thread,
        conda:
            conda_env
        shell:
            '''
            STAR {params.wasp} {params.additional} \
            --outSAMattributes CR UR {params.wasp_sam_attr} \
            --genomeDir "{input.index}" \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist None \
            --soloBarcodeMate 1 \
            --clip5pNbases 30 0 \
            --soloCBstart 1 --soloCBlen 20 \
            --soloUMIstart 21 --soloUMIlen 10 \
            --soloUMIdedup Exact --soloUMIfiltering MultiGeneUMI_All \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --soloFeatures Gene GeneFull SJ \
            --outSAMattrRGline "ID:{params.id}" "SM:{params.id}" \
            --outFileNamePrefix "{output.s1}/" \
            --readFilesIn {params.r1} {params.r2} &> "{log}"
            samtools index -@ {params.nthreads} "{params.b1}" &>> "{log}"
            '''

    rule starsolo_align_post_splitting_2:
        '''STARsolo to paired-end-align FASTQ files and generate count matrices and BAM files
           We ignore the first 30-bp (20 bp barcode + 10 bp UMI) of R1 via clip5pNbases
           This rule is run on FASTQ files that have been split by allele
        '''
        input:
            r1_split2=get_mapping_input_r1_split2,
            r2_split2=get_mapping_input_r2_split2,
            index=get_star_index,
            vcf=get_vcf_file
        output:
            s2=directory(os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2))
        threads:
            16
        log:
            os.path.join(DIR_LOGS, "{sample}.post_split2.starsolo." + align_id + ".log")
        params:
            wasp="--waspOutputMode SAMtag --varVCFfile " + vcf_file_input_param if use_vcf_wasp_workflow else "",
            wasp_sam_attr="NH HI AS nM NM NM MD vA vG vW" if use_vcf_wasp_workflow else "",
            id="{sample}",
            additional=star_additional_params,
            b1=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_, "Aligned.sortedByCoord.out.bam"),
            b2=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2, "Aligned.sortedByCoord.out.bam"),
            r1=lambda wildcards, input: ",".join([f'"{x}"' for x in input.r1_split2]),
            r2=lambda wildcards, input: ",".join([f'"{x}"' for x in input.r2_split2]),
            nthreads=lambda wildcards, input, threads: threads - 1  # Since samtools uses an additional thread as the main thread,
        conda:
            conda_env
        shell:
            '''
            STAR {params.wasp} {params.additional} \
            --outSAMattributes CR UR {params.wasp_sam_attr} \
            --genomeDir "{input.index}" \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist None \
            --soloBarcodeMate 1 \
            --clip5pNbases 30 0 \
            --soloCBstart 1 --soloCBlen 20 \
            --soloUMIstart 21 --soloUMIlen 10 \
            --soloUMIdedup Exact --soloUMIfiltering MultiGeneUMI_All \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterType BySJout \
            --soloFeatures Gene GeneFull SJ \
            --outSAMattrRGline "ID:{params.id}" "SM:{params.id}" \
            --outFileNamePrefix "{output.s2}/" \
            --readFilesIn {params.r1} {params.r2} &> "{log}"
            samtools index -@ {params.nthreads} "{params.b2}" &>> "{log}"
            '''

    rule generate_anndata_starsolo_post_splitting:
        '''Create allele-specific anndata object from STARsolo matrices
        '''
        input:
            dir_s1=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_),
            dir_s2=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2),
            mapgene=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "Gene", "raw", "barcodes.ids.txt"),
            mapgenefull=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, "Solo.out", "GeneFull", "raw", "barcodes.ids.txt")
        output:
            gene_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s1_, "Gene", adata_filename),
            genefull_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s1_, "GeneFull", adata_filename),
            gene_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s2, "Gene", adata_filename),
            genefull_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s2, "GeneFull", adata_filename)
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "{sample}.anndata.post_split.starsolo." + align_id + ".log")
        params:
            dir_s1=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_, "Solo.out", "Gene", "raw"),
            out_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s1_, "Gene"),
            dirfull_s1=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s1_, "Solo.out", "GeneFull", "raw"),
            outfull_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s1_, "GeneFull"),
            dir_s2=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2, "Solo.out", "Gene", "raw"),
            out_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s2, "Gene"),
            dirfull_s2=os.path.join(DIR_WORKUP, "results", "output_{sample}_starsolo_" + align_id, vcf_s2, "Solo.out", "GeneFull", "raw"),
            outfull_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s2, "GeneFull")
        conda:
            conda_env
        shell:
            '''
            python "{prep_anndata}" STAR "{params.dir_s1}" "{input.mapgene}" "{params.out_s1}" "" "{select_tags}" &> "{log}"
            python "{prep_anndata}" STAR "{params.dirfull_s1}" "{input.mapgenefull}" "{params.outfull_s1}" "" "{select_tags}" &>> "{log}"
            python "{prep_anndata}" STAR "{params.dir_s2}" "{input.mapgene}" "{params.out_s2}" "" "{select_tags}" &>> "{log}"
            python "{prep_anndata}" STAR "{params.dirfull_s2}" "{input.mapgenefull}" "{params.outfull_s2}" "" "{select_tags}" &>> "{log}"
            '''

    if merge_anndata_files:
        rule merge_anndata_post_splitting:
            '''Create a merged allele-specific anndata file
            '''
            input:
                gene_s1=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s1_, "Gene", adata_filename), sample=ALL_SAMPLES),
                genefull_s1=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s1_, "GeneFull", adata_filename), sample=ALL_SAMPLES),
                gene_s2=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s2, "Gene", adata_filename), sample=ALL_SAMPLES),
                genefull_s2=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "starsolo_" + align_id, vcf_s2, "GeneFull", adata_filename), sample=ALL_SAMPLES)
            output:
                gene_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "Gene", adata_filename),
                genefull_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "GeneFull", adata_filename),
                gene_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "Gene", adata_filename),
                genefull_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "GeneFull", adata_filename),
                info_gene_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "Gene", "sample_order.txt"),
                info_genefull_s1=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s1_, "GeneFull", "sample_order.txt"),
                info_gene_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "Gene", "sample_order.txt"),
                info_genefull_s2=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "starsolo_" + align_id, vcf_s2, "GeneFull", "sample_order.txt")
            log:
                os.path.join(DIR_LOGS, "anndata.merged.post_split." + align_id + ".log")
            threads:
                1
            conda:
                conda_env
            shell:
                '''
                python "{merge_anndata}" "{output.gene_s1}" {input.gene_s1} &> "{log}"
                python "{merge_anndata}" "{output.genefull_s1}" {input.genefull_s1} &>> "{log}"
                python "{merge_anndata}" "{output.gene_s2}" {input.gene_s2} &>> "{log}"
                python "{merge_anndata}" "{output.genefull_s2}" {input.genefull_s2} &>> "{log}"
                (echo "{input.gene_s1}" > "{output.info_gene_s1}") &>> "{log}"
                (echo "{input.genefull_s1}" > "{output.info_genefull_s1}") &>> "{log}"
                (echo "{input.gene_s2}" > "{output.info_gene_s2}") &>> "{log}"
                (echo "{input.genefull_s2}" > "{output.info_genefull_s2}") &>> "{log}"
                '''


rule generate_index_star:
    '''Create a STAR genome index from FASTA and GTF
    '''
    input:
        fasta=genome_fasta,
        gtf=genome_gtf
    output:
        index=directory(star_index_dir_new),
        fasta_unzipped=temp(os.path.join(DIR_WORKUP, "tmp", "temp.fasta")),
        gtf_unzipped=temp(os.path.join(DIR_WORKUP, "tmp", "temp.gtf"))
    threads:
        14
    log:
        os.path.join(DIR_LOGS, "star.index" + ".log")
    conda:
        conda_env
    shell:
        '''
        zcat -f "{input.fasta}" > "{output.fasta_unzipped}"
        zcat -f "{input.gtf}" > "{output.gtf_unzipped}"
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir "{output.index}" --sjdbOverhang 100 \
        --genomeFastaFiles "{output.fasta_unzipped}" --sjdbGTFfile "{output.gtf_unzipped}" &> "{log}"
        '''



rule kallisto_align:                                                            
    '''kallisto to map paired-end-align FASTQ files and generate nascent/mature/ambiguous count matrices
    '''                                                                         
    input:
        r1=get_mapping_input_r1,
        r2=get_mapping_input_r2,
        i=get_kallisto_index,     
        g=get_kallisto_t2g,
        c1=get_kallisto_c1,
        c2=get_kallisto_c2,
        map=os.path.join(DIR_WORKUP, "assigned", "{sample}_mapping.txt")
    output:
        directory(os.path.join(DIR_WORKUP, "results", "output_{sample}_kallisto_" + align_id, "counts_unfiltered"))
    threads:
        4*thread_factor
    log:
        os.path.join(DIR_LOGS, "{sample}.kallisto." + align_id + ".log")
    params:                               
        id="{sample}",                                               
        out=os.path.join(DIR_WORKUP, "results", "output_{sample}_kallisto_" + align_id),
        barcodes=os.path.join(DIR_WORKUP, "results", "output_{sample}_kallisto_" + align_id, "counts_unfiltered", "cells_x_genes.barcodes.txt"),
        barcodeIDs=os.path.join(DIR_WORKUP, "results", "output_{sample}_kallisto_" + align_id, "counts_unfiltered", "cells_x_genes.barcodes.ids.txt"),
        reads=lambda wildcards, input: " ".join(
            [f'"{file}"' for pair in zip(input.r1, input.r2) for file in pair]
        )
    conda:
        conda_env                     
    shell:      
        '''
        kb count --workflow=nac --sum=total -x 0,0,20:0,20,30:0,30,0,1,0,0 \
        -i "{input.i}" \
        -g "{input.g}" \
        -c1 "{input.c1}" \
        -c2 "{input.c2}" \
        -t {threads} \
        --parity=paired \
        --strand=forward \
        -w None \
        -o "{params.out}/" \
        {params.reads} &> "{log}"
        (python "{decode_barcode}" "{input.map}" "{params.barcodes}" "{params.barcodeIDs}") &>> "{log}"
        '''                                         


rule generate_index_kallisto:
    '''Create a kallisto transcriptome index from FASTA and GTF
    '''
    input:
        fasta=genome_fasta,
        gtf=genome_gtf
    output:
        i=kallisto_index_new,
        g=kallisto_t2g_new,
        c1=kallisto_c1_new,
        c2=kallisto_c2_new
    threads:
        14
    log:
        os.path.join(DIR_LOGS, "kallisto.index" + ".log")
    params:
        f1=kallisto_f1_new,
        f2=kallisto_f2_new,
        dlist=kallisto_dlist
    conda:
        conda_env
    shell:
        '''
        kb ref --workflow=nac -t {threads} -i "{output.i}" -g "{output.g}" \
        -c1 "{output.c1}" -c2 "{output.c2}" \
        -f1 "{params.f1}" -f2 "{params.f2}" \
        --d-list="{params.dlist}" \
        "{input.fasta}" "{input.gtf}" &> "{log}"
        '''

if use_features:
    rule kallisto_index_features:
        '''Use kallisto to index features from FASTA file
        '''
        input:
            features_fasta=features_fasta_file
        output:
            i=temp(os.path.join(DIR_WORKUP, "kallisto", "index.idx")),
            g=temp(os.path.join(DIR_WORKUP, "kallisto", "t2g.txt"))
        threads:
            4
        params:
            k=features_kmer_length
        log:
            os.path.join(DIR_LOGS, "kallisto.index.features.log")
        conda:
            conda_env
        shell:
            '''
            (paste \
            <(zcat -f "{input.features_fasta}"|grep "^>"|sed 's/^>//'|cut -d' ' -f1) \
            <(zcat -f "{input.features_fasta}"|grep "^>"|sed 's/^>//'|cut -d' ' -f1) > "{output.g}") &> "{log}"
            kb ref --overwrite --workflow=custom -t {threads} -i "{output.i}" \
            -k {params.k} "{input.features_fasta}" &>> "{log}"
            '''

    rule kallisto_map_features:
        '''Create a barcode-x-feature count matrix by mapping reads to features
        '''
        input:
            r1=get_mapping_input_r1,
            r2=get_mapping_input_r2,
            i=os.path.join(DIR_WORKUP, "kallisto", "index.idx"),
            g=os.path.join(DIR_WORKUP, "kallisto", "t2g.txt"),
            map=os.path.join(DIR_WORKUP, "assigned", "{sample}_mapping.txt")
        output:
            file=os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.mtx"),
            barcodeIDs=os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.barcodes.ids.txt")
        threads:
            4*thread_factor
        params:
            dir=os.path.join(DIR_WORKUP, "results", "output_{sample}_features"),
            k=features_kmer_length,
            start=20+10+features_start,
            end=-1 if features_length <= 0 else 20+10+features_start+features_length,
            barcodes=os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.barcodes.txt"),
            reads=lambda wildcards, input: " ".join([f'"{file}"' for file in input.r1])
        log:
            os.path.join(DIR_LOGS, "{sample}.kallisto.features.log")
        conda:
            conda_env
        shell:
            '''
            kb count --overwrite -t {threads} -i "{input.i}" -g "{input.g}" \
            -w None -x "0,0,20:0,20,30:0,{params.start},{params.end}" -o "{params.dir}" {params.reads} &> "{log}"
            touch "{output.file}" && touch "{params.barcodes}" &>> "{log}"
            (python "{decode_barcode}" "{input.map}" "{params.barcodes}" "{output.barcodeIDs}") &>> "{log}"
            '''

    rule generate_anndata_kallisto_features:
        '''Create anndata object from kallisto feature matrices
        '''
        input:
            mtx=os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.mtx"),
            map=os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered", "cells_x_genes.barcodes.ids.txt")
        output:
            os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "features", adata_filename)
        threads:
            1
        log:
            os.path.join(DIR_LOGS, "{sample}.anndata.kallisto.features.log")
        params:
            dir=os.path.join(DIR_WORKUP, "results", "output_{sample}_features", "counts_unfiltered"),
            out=os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "features")
        conda:
            conda_env
        shell:
            '''
            python "{prep_anndata}" kallisto "{params.dir}" "{input.map}" "{params.out}" "" "{select_tags}" &> "{log}"
            '''



if use_features and merge_anndata_files:
    rule merge_features_anndata:
        '''Create a merged anndata file from features
        '''
        input:
            feature=expand(os.path.join(DIR_WORKUP, "results", "anndatas", "{sample}_anndata", "features", adata_filename), sample=ALL_SAMPLES)
        output:
            feature=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "features", adata_filename),
            info_feature=os.path.join(DIR_WORKUP, "results", "anndatas", "merged", "features", "sample_order.txt")
        log:
            os.path.join(DIR_LOGS, "anndata.kallisto.features.merged." + align_id + ".log")
        threads:
            1
        conda:
            conda_env
        shell:
            '''
            python "{merge_anndata}" "{output.feature}" {input.feature} &> "{log}"
            (echo "{input.feature}" > "{output.info_feature}") &>> "{log}"
            '''


onsuccess:
    shell(f'python "{html_report_gen}" "{DIR_WORKUP}" {align_id} > {DIR_LOGS}/final_script.log 2>&1')


