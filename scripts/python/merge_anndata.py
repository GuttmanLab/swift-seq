# Script: merge_anndata.py
# Author: Delaney K. Sullivan
# Usage: python merge_anndata.py output_file.h5ad anndata_file1.h5ad anndata_file2.h5ad ... [barcodes_to_keep (optional)]

import os
import os.path
import sys

# Check command-line arguments

if len(sys.argv) < 3:
    print("")
    print("Usage: python merge_anndata.py output_file.h5ad anndata_file1.h5ad anndata_file2.h5ad ... [barcodes_to_keep (optional)]")
    print("")
    print("(barcodes_to_keep, an optional argument, can be something like RTBC1,RTBC3 which means only barcode IDs ending in the [RTBC1] or the [RTBC3] tags will be retained)")
    print("")
    sys.exit(1)


# Some basic settings

out_fname = sys.argv[1] # Output filename
arg_filter = None # What filter should be applied to the barcode ID
anndata_files = sys.argv[2:]
final_arg = sys.argv[len(sys.argv)-1]
if not os.path.isfile(final_arg):
    arg_filter = final_arg # Set the filter
    anndata_files = anndata_files[:-1] # Remove last element, since that's the filter, not a file


# Import libraries

import numpy as np
import scipy
import anndata
import random
import pandas as pd
import warnings

def nd(arr):
    return np.asarray(arr).reshape(-1)

def load_and_concatenate_anndata(files):
    # Load the first Anndata object to initialize the concatenation process
    i = 0
    adata_combined = None
    init = False
    for filename in files[0:]:
        adata_combined = anndata.read_h5ad(files[i])
        i = i + 1
        rows, columns = adata_combined.shape
        if rows != 0 and columns != 0:
            init = True
            break

    if len(files) == 1 or not init:
        return adata_combined

    adata_combined.obs.index = [x + f'-{i}' for x in adata_combined.obs.index]
    # Iterate over remaining Anndata files and concatenate them one by one
    start = i
    for filename in files[start:]:
        i = i+1
        adata = anndata.read_h5ad(filename)
        rows, columns = adata.shape
        if rows == 0 or columns == 0:
            continue
        adata.obs.index = [x + f'-{i}' for x in adata.obs.index]
        adata_combined = anndata.concat([adata_combined, adata], axis=0, join='outer', merge='same')

    return adata_combined


warnings.filterwarnings('ignore')

# Merge the anndata objects

data = load_and_concatenate_anndata(anndata_files)

# Check if anndata looks good

rows, columns = data.shape
if rows == 0 or columns == 0:
    print("Note: Empty AnnData")
    data.write_h5ad(out_fname)
    sys.exit(0)

# Barcode filtering

if arg_filter:
    filter_list = ['[' + item + ']' for item in arg_filter.split(',')]
    data = data[data.obs['id'].str.endswith(tuple(filter_list))]

data.obs["cell_counts"] = data.X.sum(axis=1)
data.var["gene_counts"] = nd(data.X.sum(axis=0))
data.obs["n_genes"] = nd((data.X>0).sum(axis=1))
data.var["n_cells"] = nd((data.X>0).sum(axis=0))
mito_genes = data.var["gene_name"].str.lower().str.startswith('mt-')
data.obs["percent_mito"] = data[:,mito_genes].X.sum(axis=1)/data.X.sum(axis=1)*100

print(data)

# Save anndata object to file

data.write_h5ad(out_fname)
print(f"Concatenated Anndata object saved to {out_fname}")

