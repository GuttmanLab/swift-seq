import sys

# Some basic settings

program = sys.argv[1] # STAR or kallisto or kallisto_full
counts_dir = sys.argv[2] # "./out_dir/workup/results/output_PBMC_starsolo_hg38_masked/Solo.out/Gene/raw"
barcode_mapping_file = sys.argv[3] # "./out_dir/workup/assigned/PBMC_mapping.txt"
output_folder = sys.argv[4] # "./out_dir/workup/results/output_PBMC_starsolo_hg38_masked"
output_prefix = sys.argv[5]
filter_for_tags = sys.argv[6] # "RTBC1,RTBC2,RTBC3,RTBC4"

# For STAR, the counts_dir should be something like: output/Solo.out/Gene/raw
# (Use output/Solo.out/GeneFull/raw instead to also include introns in the count quantifications)

# Kallisto options (if using kallisto cells_x_genes mtx files)
## The output directory specified should be something like output/counts_unfiltered

use_kallisto = False # Set to true if using kallisto matrices
kallisto_geneFull = False # Set to true to also count introns

if program.lower() == "kallisto":
    use_kallisto = True
if program.lower() == "kallisto_full":
    use_kallisto = True
    kallisto_geneFull = True

# Import libraries

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

import scipy
import scipy.io
import anndata
import random
import pandas as pd
import warnings
import matplotlib.patches as mpatches

def nd(arr):
    return np.asarray(arr).reshape(-1)

# Violin plot functions

def vplot(y, ax):
    parts = ax.violinplot(
        y, showmeans=False, showmedians=False,
        showextrema=False)

    mean = y.mean()
    ax.scatter(1, mean, zorder=10, color="white")
    
    x = np.random.normal(1, 0.04, size=len(y))
    ax.scatter(x, y, color="k", s=1)
    
    for pc in parts['bodies']:
        pc.set_facecolor('#D43F3A')
        pc.set_edgecolor('black')
        pc.set_alpha(1)
    
    ax.set_xticks([1])
    ax.set_xticklabels([""])
    return ax

warnings.filterwarnings('ignore')

fsize=20
plt.rcParams.update({'font.size': fsize})

# STAR input files

matrix_file="matrix.mtx"
barcodes_file="barcodes.tsv"
genes_file="features.tsv"
genes_names_file="features.tsv"

# Some kallisto-specific options

if use_kallisto:
    genes_file="cells_x_genes.genes.txt"
    genes_names_file="cells_x_genes.genes.names.txt"
    barcodes_file="cells_x_genes.barcodes.txt"
    if kallisto_geneFull:
        matrix_file="cells_x_genes.total.mtx"
    else:
        matrix_file="cells_x_genes.cell.mtx"
        if not os.path.exists(os.path.join(counts_dir, matrix_file)):
            matrix_file="cells_x_genes.mtx"

# Set up paths

matrix_path=os.path.join(counts_dir, matrix_file)
barcodes_path=os.path.join(counts_dir, barcodes_file)
genes_path=os.path.join(counts_dir, genes_file)
genes_names_path=os.path.join(counts_dir, genes_names_file)

# If matrix file exists but is empty, just write an empty AnnData

if os.path.exists(matrix_path) and os.stat(matrix_path).st_size == 0:
    empty_adata = anndata.AnnData()
    empty_adata.write(os.path.join(output_folder, output_prefix + "adata.h5ad"))
    sys.exit(0)

# Read in matrix, genes, and barcodes

with open(genes_path, 'r') as f:
    genes = [line.strip().split()[0] for line in f]

with open(genes_names_path, 'r') as f:
    if genes_path == genes_names_path:
        genes_names = [line.strip().split()[1] for line in f]
    else:
        genes_names = [line.strip().split()[0] for line in f]

with open(barcodes_path, 'r') as f:
    barcodes = [line.strip() for line in f]

data = anndata.AnnData(scipy.io.mmread(matrix_path).tocsr()).T  # float32 matrix
if use_kallisto:
    data = data.T
data.obs_names = barcodes
data.var_names = genes
data.var['gene_id'] = genes
data.var['gene_name'] = genes_names

# Barcode mapping and filtering

mapping_df = pd.read_csv(barcode_mapping_file, sep='\t', header=None, names=['obs_index', 'new_data'])
mapping_df['new_data'] = mapping_df['new_data'].apply(lambda x: '[' + x.replace(',', '][') + ']')
mapping = pd.Series(mapping_df.new_data.values, index=mapping_df.obs_index).to_dict()
data.obs['id'] = data.obs.index.map(mapping)
filter_for_tags_final = ['[' + item + ']' for item in filter_for_tags.split(',')]
if filter_for_tags.strip():
    filter_criteria = pd.Series([False] * len(mapping_df), index=mapping_df.index)
    for tag in filter_for_tags_final:
        filter_criteria |= mapping_df['new_data'].str.strip().str.contains(tag, regex=False)
    filtered_df = mapping_df[filter_criteria]
    barcodes_to_keep = filtered_df['obs_index'].unique().tolist()
    data = data[data.obs_names.isin(barcodes_to_keep), :]

data.obs["cell_counts"] = data.X.sum(axis=1)
data.var["gene_counts"] = nd(data.X.sum(axis=0))
data.obs["n_genes"] = nd((data.X>0).sum(axis=1))
data.var["n_cells"] = nd((data.X>0).sum(axis=0))
mito_genes = data.var["gene_name"].str.lower().str.startswith('mt-')
data.obs["percent_mito"] = data[:,mito_genes].X.sum(axis=1)/data.X.sum(axis=1)*100
data.var_names = list(data.var.gene_name)
data.var_names_make_unique()




print(data)


# Save anndata object to file

data.write_h5ad(os.path.join(output_folder, output_prefix + "adata.h5ad"))
rows, columns = data.shape
if rows == 0 or columns == 0:
    print("Note: The AnnData object has zero rows or columns.")
    sys.exit(0)

# Make QC plots

fig, ax = plt.subplots(figsize=(5, 5))
x = nd(data.X.sum(axis=1))
y = nd(np.sum(data.X>0, axis=1))
ax.scatter(x, y, color="green", alpha=0.25)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
ax.set_xlim(0)
ax.set_ylim(0)
plt.savefig(os.path.join(output_folder, output_prefix + "counts_scatter.png"))


knee = np.sort(nd(data.X.sum(axis=1)))[::-1]

fig, ax = plt.subplots(figsize=(5, 5))

x = knee
y = range(len(knee))

ax.loglog(x, y, linewidth=5, color="g")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.savefig(os.path.join(output_folder, output_prefix + "knee.png"))


fig, ax = plt.subplots(figsize=(5*3,5), ncols=3)
x1 = data.obs["n_genes"]
x2 = nd(data.X.sum(axis=1))
x3 = data.obs["percent_mito"]
vplot(x1, ax[0])
vplot(x2, ax[1])
vplot(x3, ax[2])
ax[0].set_ylabel("Genes detected")
ax[1].set_ylabel("UMI counts")
ax[2].set_ylabel("Percent Mito")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, output_prefix + "counts_violin.png"))


