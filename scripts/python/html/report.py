#!/usr/bin/env python3

import os
import sys
import glob
import re
import base64
import json
import anndata
import scanpy as sc
import warnings
import csv
import pandas as pd
from scipy.io import mmread

warnings.filterwarnings('ignore')


from knee_plots import create_knee_plot_html
from knee_plots import create_knee_plot_html_from_anndata
from knee_plots import create_sequencing_saturation_plot
from knee_plots import create_mouse_human_scatter
from knee_plots import create_nascent_mature_ambiguous_bars
from knee_plots import create_suffix_umi_comparison_box
from knee_plots import _PLOTLY_JS

import scipy.sparse as sp
import numpy as np


###################################################################
# FASTQC toggles
###################################################################

def read_fastqc_html(path):
    if not os.path.exists(path):
        return None
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        return f.read()

def generate_fastqc_iframe_toggles(sample_name, base_dir):
    r1_path= os.path.join(base_dir,"qc",f"{sample_name}_R1_fastqc.html")
    r2_path= os.path.join(base_dir,"qc",f"{sample_name}_R2_fastqc.html")

    r1_html= read_fastqc_html(r1_path)
    r2_html= read_fastqc_html(r2_path)

    lines=[]
    lines.append("<div class='fastqc-toggles' style='margin:1em 0;'>")

    if r1_html:
        lines.append(f"""
<div style="margin-bottom:0.5em;">
  <button onclick="toggleFastQC('{sample_name}_R1_div')">Toggle R1 FastQC</button>
  <div id="{sample_name}_R1_div" style="display:none; margin-top:0.5em;">
    <iframe width="100%" height="800px" srcdoc="{r1_html.replace('"','&quot;').replace("'",'&#39;')}"></iframe>
  </div>
</div>
""")

    if r2_html:
        lines.append(f"""
<div style="margin-bottom:0.5em;">
  <button onclick="toggleFastQC('{sample_name}_R2_div')">Toggle R2 FastQC</button>
  <div id="{sample_name}_R2_div" style="display:none; margin-top:0.5em;">
    <iframe width="100%" height="800px" srcdoc="{r2_html.replace('"','&quot;').replace("'",'&#39;')}"></iframe>
  </div>
</div>
""")

    lines.append("</div>")
    return "\n".join(lines)


###################################################################
# Parsing logs for alignment, etc.
###################################################################

def parse_summary_csv(base_dir, sample, run_name, which):
    """
    which should be "GeneFull" or "Gene".
    Returns a dict mapping keys → float.
    """
    path = os.path.join(
        base_dir, "results",
        f"output_{sample}_starsolo_{run_name}",
        "Solo.out", which, "Summary.csv"
    )
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 2: 
                continue
            key, val = row
            try:
                out[key] = float(val)
            except ValueError:
                pass
    return out

def parse_log_file(log_path):
    metrics = {
        "input_reads": 0,
        "avg_input_read_length": 0.0,
        "unique_pct": 0.0,
        "multimapped_pct": 0.0,
        "multimapped_too_many_loci_pct": 0.0,
        "unmapped_too_many_mismatches_pct": 0.0,
        "unmapped_too_short_pct": 0.0,
        "unmapped_other_pct": 0.0,
        "unmapped_other_count": 0,
        "chimeric_pct": 0.0,
        "chimeric_count": 0,
        "avg_mapped_length": 0.0,
        "num_splices_total": 0,
        "num_splices_annotated": 0
    }
    if not os.path.exists(log_path):
        return metrics
    
    with open(log_path,"r") as f:
        for line in f:
            ln=line.strip()
            if ln.startswith("Number of input reads |"):
                parts= ln.split("|")
                metrics["input_reads"]= int(parts[1].strip())
            elif ln.startswith("Average input read length |"):
                parts= ln.split("|")
                metrics["avg_input_read_length"]= float(parts[1].strip())
            elif ln.startswith("Uniquely mapped reads % |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["unique_pct"]= float(pct_str)
            elif ln.startswith("% of reads mapped to multiple loci |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["multimapped_pct"]= float(pct_str)
            elif ln.startswith("% of reads mapped to too many loci |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["multimapped_too_many_loci_pct"]= float(pct_str)
            elif ln.startswith("% of reads unmapped: too many mismatches |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["unmapped_too_many_mismatches_pct"]= float(pct_str)
            elif ln.startswith("% of reads unmapped: too short |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["unmapped_too_short_pct"]= float(pct_str)
            elif ln.startswith("% of reads unmapped: other |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["unmapped_other_pct"]= float(pct_str)
            elif ln.startswith("Number of reads unmapped: other |"):
                parts= ln.split("|")
                metrics["unmapped_other_count"]= int(parts[1].strip())
            elif ln.startswith("Number of chimeric reads |"):
                parts= ln.split("|")
                metrics["chimeric_count"]= int(parts[1].strip())
            elif ln.startswith("% of chimeric reads |"):
                parts= ln.split("|")
                pct_str= parts[1].strip().replace("%","")
                metrics["chimeric_pct"]= float(pct_str)
            elif ln.startswith("Average mapped length |"):
                parts= ln.split("|")
                metrics["avg_mapped_length"]= float(parts[1].strip())
            elif ln.startswith("Number of splices: Total |"):
                parts= ln.split("|")
                metrics["num_splices_total"]= int(parts[1].strip())
            elif ln.startswith("Number of splices: Annotated (sjdb) |"):
                parts= ln.split("|")
                metrics["num_splices_annotated"]= int(parts[1].strip())
    return metrics

def parse_trimmed_qc_files(base_dir, sample_name):
    import re
    pat= re.compile(r"^Pairs written \(passing filters\):\s+([\d,]+)")
    tot=0
    pattern= os.path.join(base_dir,"trimmed",f"{sample_name}_*.trimmed.qc.txt")
    fls= glob.glob(pattern)
    for fpath in fls:
        if os.path.isfile(fpath):
            with open(fpath,"r") as f:
                for line in f:
                    m= pat.search(line.strip())
                    if m:
                        vs= m.group(1).replace(",","")
                        tot+= int(vs)
    return tot

def parse_kallisto_info(base_dir, sample_name, run_name):
    """
    If out_dir/workup/results/output_<sample>_kallisto_<run_name>/run_info.json exists,
    return a dict with:
      - n_processed
      - n_pseudoaligned
      - n_unique
      - p_pseudoaligned (percent)
      - p_unique      (percent)
    Otherwise return None.
    """
    path = os.path.join(
        base_dir, "results",
        f"output_{sample_name}_kallisto_{run_name}",
        "run_info.json"
    )
    if not os.path.exists(path):
        return None

    with open(path) as f:
        info = json.load(f)

    # absolute counts
    n_proc   = info.get("n_processed",    0)
    n_pseudo = info.get("n_pseudoaligned", 0)
    n_unique = info.get("n_unique",        0)

    # computed percentages
    p_pseudo = (n_pseudo / n_proc * 100.0) if n_proc else 0.0
    p_unique = (n_unique / n_proc * 100.0) if n_proc else 0.0

    call_str = info.get("call", "")

    return {
        "n_processed":      n_proc,
        "n_pseudoaligned":  n_pseudo,
        "n_unique":         n_unique,
        "p_pseudoaligned":  p_pseudo,
        "p_unique":         p_unique,
        "call":             call_str
    }

def parse_cell_reads_stats(base_dir, sample_name, run_name):
    """
    Parse $base_dir/results/output_<sample>_starsolo_<run_name>/Solo.out/GeneFull/CellReads.stats
    to get total reads (cbMatch) and unique UMIs (nUMIunique) per cell.
    Return a list of (reads, umis) so we can compute mean reads / cell and saturation, etc.
    If the file doesn't exist, return an empty list.
    """
    stats_path = os.path.join(
        base_dir, "results",
        f"output_{sample_name}_starsolo_{run_name}",
        "Solo.out", "GeneFull", "CellReads.stats"
    )
    if not os.path.exists(stats_path):
        return []

    rows = []
    with open(stats_path, "r") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            # skip header
            if header is None:
                header = line.split()
                continue
            parts = line.split()
            cbMatch_index = header.index("cbMatch")
            nUMIunique_index = header.index("nUMIunique")
            reads_val = float(parts[cbMatch_index])
            umis_val = float(parts[nUMIunique_index])
            rows.append((reads_val, umis_val))

    return rows

def parse_ligation_file(ligation_path):
    result= {"barcode_dict":{}, "total_reads_post_trimming":0}
    if not os.path.exists(ligation_path):
        return result
    lines_before_blank=[]
    with open(ligation_path,"r") as f:
        for line in f:
            ln=line.strip()
            if ln=="":
                break
            lines_before_blank.append(ln)
    import re
    pat= re.compile(r"^(\d+)\s*\(.*?\)\s*reads found with\s+(\d+)\s+barcodes?")
    for l in lines_before_blank:
        mm= pat.search(l)
        if mm:
            cval= int(mm.group(1))
            bc= int(mm.group(2))
            result["barcode_dict"][bc]= result["barcode_dict"].get(bc,0)+ cval
    total= sum(result["barcode_dict"].values())
    result["total_reads_post_trimming"]= total
    return result

def merge_ligation_dicts(lst):
    merged_bc={}
    sum_all=0
    for ld in lst:
        for bc_c, val in ld["barcode_dict"].items():
            merged_bc[bc_c]= merged_bc.get(bc_c,0)+ val
        sum_all+= ld["total_reads_post_trimming"]
    return {"barcode_dict": merged_bc, "total_reads_post_trimming": sum_all}

def parse_exclusion_file(exclusion_path):
    result={"excluded_count":0,"excluded_total":0,"excluded_pct":0.0}
    if not os.path.exists(exclusion_path):
        return result
    lines=[]
    with open(exclusion_path,"r") as f:
        for ln in f:
            lnn= ln.strip()
            if lnn:
                lines.append(lnn)
    if len(lines)<2:
        return result
    import re
    pat1= re.compile(r"^(\d+)\s+aligned to exclusion index out of\s+(\d+)\s+reads")
    pat2= re.compile(r"^Percent excluded:\s+(\d+(\.\d+)?)%")
    m1= pat1.search(lines[0])
    m2= pat2.search(lines[1])
    if m1:
        result["excluded_count"]= int(m1.group(1))
        result["excluded_total"]= int(m1.group(2))
    if m2:
        result["excluded_pct"]= float(m2.group(1))
    return result

def merge_exclusion_dicts(lst):
    te=0
    tr=0
    for e in lst:
        te+= e["excluded_count"]
        tr+= e["excluded_total"]
    merged={"excluded_count":te,"excluded_total":tr,"excluded_pct":0.0}
    if tr>0:
        merged["excluded_pct"]=(te/tr)*100.0
    return merged

def parse_fastqc_total_reads(fastqc_html_path):
    """
    Fixing the regex to match the actual text.
    """
    if not os.path.exists(fastqc_html_path):
        return 0
    import re
    # The correct pattern (no double backslashes):
    pat= re.compile(r"<td>\s*Total\s+Sequences\s*</td>\s*<td>\s*(\d+)\s*</td>")
    with open(fastqc_html_path,"r",encoding="utf-8",errors="ignore") as f:
        txt= f.read()
    mm= pat.search(txt)
    if mm:
        return int(mm.group(1))
    return 0

def read_anndata_objects(base_dir, sample_name, run_name, return_merge=False):
    """
    Loads two AnnData objects for the sample:
      - GeneFull (with introns)
      - Gene     (excluding introns)
    Returns (adata_gene_full, adata_gene).
    If a file doesn't exist, you can return None for that.
    """
    import os
    import anndata

    # Build paths
    sname = "merged" if return_merge else f"{sample_name}_anndata"
    gf_path = os.path.join(
        base_dir, "results", "anndatas",
        sname,
        f"starsolo_{run_name}", "GeneFull", "adata.h5ad"
    )
    gene_path = os.path.join(
        base_dir, "results", "anndatas",
        sname,
        f"starsolo_{run_name}", "Gene", "adata.h5ad"
    )

    adata_gf = None
    adata_g  = None

    if os.path.exists(gf_path):
        adata_gf = anndata.read_h5ad(gf_path)
        adata_gf.X = adata_gf.X.astype(np.int32)
        adata_gf.obs_names = [f"{bar}_{sample_name}" for bar in adata_gf.obs_names]
    if os.path.exists(gene_path):
        adata_g  = anndata.read_h5ad(gene_path)
        adata_g.X = adata_g.X.astype(np.int32)
        adata_g.obs_names = [f"{bar}_{sample_name}" for bar in adata_g.obs_names]

    return adata_gf, adata_g


def load_kallisto_anndata(base_dir: str, sample_name: str, run_name: str):
    """
    Load Kallisto count matrices into two AnnData objects:
      - GeneFull (including introns) from counts_unfiltered/cells_x_genes.total.mtx
      - Gene     (excluding introns) from counts_unfiltered/cells_x_genes.cell.mtx

    Returns:
        adata_total, adata_cell
        (either or both may be None if the corresponding file doesn't exist)
    """
    prefix = os.path.join(
        base_dir, "results",
        f"output_{sample_name}_kallisto_{run_name}",
        "counts_unfiltered"
    )
    total_path = os.path.join(prefix, "cells_x_genes.total.mtx")
    cell_path  = os.path.join(prefix, "cells_x_genes.cell.mtx")

    adata_total = None
    adata_cell  = None

    # GeneFull
    if os.path.exists(total_path):
        mat = mmread(total_path).tocsr()
        adata_total = anndata.AnnData(X=mat)
        # optional: give meaningful obs/var names if you have them
        adata_total.obs_names = [f"cell{i}" for i in range(mat.shape[0])]
        adata_total.var_names = [f"gene{i}" for i in range(mat.shape[1])]

    # Gene
    if os.path.exists(cell_path):
        mat = mmread(cell_path).tocsr()
        adata_cell = anndata.AnnData(X=mat)
        adata_cell.obs_names = [f"cell{i}" for i in range(mat.shape[0])]
        adata_cell.var_names = [f"gene{i}" for i in range(mat.shape[1])]

    if not adata_total and not adata_cell:
        return None, None

    # Update the anndata
    with open(os.path.join(prefix, "cells_x_genes.barcodes.txt"), 'r') as f:
        barcodes = [line.strip() for line in f]
    with open(os.path.join(prefix, "cells_x_genes.genes.txt"), 'r') as f:
        genes = [line.strip() for line in f]
    with open(os.path.join(prefix, "cells_x_genes.genes.names.txt"), 'r') as f:
        genes_names = [line.strip() for line in f]

    adata_total.obs_names = barcodes
    adata_total.var_names = genes
    adata_total.var['gene_id'] = genes
    adata_total.var['gene_name'] = genes_names
    adata_cell.obs_names = barcodes
    adata_cell.var_names = genes
    adata_cell.var['gene_id'] = genes
    adata_cell.var['gene_name'] = genes_names


    # Barcode mapping

    barcode_mapping_file = os.path.join(prefix, "cells_x_genes.barcodes.ids.txt")
    mapping_df = pd.read_csv(barcode_mapping_file, sep='\t', header=None, names=['obs_index', 'new_data'])
    mapping_df['new_data'] = mapping_df['new_data'].apply(lambda x: '[' + x.replace(',', '][') + ']')
    mapping = pd.Series(mapping_df.new_data.values, index=mapping_df.obs_index).to_dict()

    adata_total.obs['id'] = adata_total.obs.index.map(mapping)
    adata_cell.obs['id'] = adata_cell.obs.index.map(mapping)


    return adata_total, adata_cell



def create_merged_anndata(all_samples_anndata, all_samples):
    """
    Efficiently merges a list of AnnData objects from a dict.
    
    Parameters:
    - all_samples_anndata: dict mapping sample names to AnnData objects (already suffixed).
    - all_samples: list of sample names to include in the merge.
    
    Returns:
    - A merged AnnData object, or None if no valid samples exist.
    """
    adata_list = [all_samples_anndata[s] for s in all_samples if s in all_samples_anndata]
    if not adata_list:
        return None
    if len(adata_list) == 1:
        return adata_list[0]  # no .copy(), assuming caller won't mutate in place
    if not adata_list[0]:
        return None
    return adata_list[0].concatenate(*adata_list[1:], join="outer", index_unique=None)

def bullet_pct(label, pct, input_reads, raw_count=None, threshold=0.1):
    if pct< threshold:
        if raw_count is None:
            raw_count= round((pct/100.0)* input_reads)
        return f"{label}: {pct:.2f}% ({raw_count:,} reads)"
    else:
        return f"{label}: {pct:.2f}%"

def parse_mtx_pseudobulk(mtx_path):
    """
    Reads a MatrixMarket .mtx of shape (genes × cells) and returns:
      - total_counts: sum of all values
      - per_cell_counts: dict {cell_index: sum of that column}
    Assumes 1‑based indices in the file.
    """
    total = 0
    per_cell = {}
    with open(mtx_path) as f:
        # skip comments
        for line in f:
            if line.startswith('%'):
                continue
            parts = line.split()
            if len(parts)==3:
                # header line: num_rows, num_cols, num_entries
                n_rows, n_cols, n_entries = map(int, parts)
                break
        # now read entries
        for _ in range(n_entries):
            row_str, col_str, val_str = f.readline().split()
            col = int(col_str) - 1  # zero‑based
            val = float(val_str)
            total += val
            per_cell[col] = per_cell.get(col, 0.0) + val
    return total, per_cell


###################################################################
# Boxes
###################################################################

def generate_read_processing_html(
    fastqc_reads,
    trimmed_pairs,
    ligation_data,
    exclusion_data,
    final_num_reads
):
    rows=[]
    bc_dict= ligation_data["barcode_dict"]
    AAA= ligation_data["total_reads_post_trimming"]
    has_exclusion= (exclusion_data["excluded_total"]>0)
    any_data= bool(fastqc_reads>0 or trimmed_pairs>0 or bc_dict or has_exclusion or final_num_reads>0)

    if not any_data:
        return ""

    if fastqc_reads>0:
        rows.append(f"<tr><td>Total reads:</td><td>{fastqc_reads:,}</td></tr>")
    if trimmed_pairs>0:
        rows.append(f"<tr><td>Number of reads post-trimming:</td><td>{trimmed_pairs:,}</td></tr>")

    if bc_dict:
        sorted_bc= sorted(bc_dict.keys())
        if AAA>0:
            for bc in sorted_bc:
                cval= bc_dict[bc]
                pct_bc= (cval/AAA)*100.0
                rows.append(f"<tr><td>Reads with {bc} barcodes (%):</td><td>{pct_bc:.2f}%</td></tr>")
            max_bc= max(sorted_bc)
            rows.append(f"<tr><td>Reads with all {max_bc} barcodes:</td><td>{bc_dict[max_bc]:,}</td></tr>")
        else:
            for bc in sorted_bc:
                cval= bc_dict[bc]
                rows.append(f"<tr><td>Reads with {bc} barcodes:</td><td>{cval:,}</td></tr>")

    ddd_note=""
    if AAA>0:
        ddd= trimmed_pairs - AAA
        if ddd>0:
            ddd_note= f"<p><em>(Note: An additional {ddd:,} reads were filtered out during barcode identification)</em></p>"

    footnote_html=""
    if has_exclusion:
        ec= exclusion_data["excluded_count"]
        et= exclusion_data["excluded_total"]
        ep= exclusion_data["excluded_pct"]
        e_str= f"{ep:.2f}% <i>({ec:,} out of {et:,} reads)</i>"
        rows.append(f"<tr><td>Reads excluded*:</td><td>{e_str}</td></tr>")
        footnote_html= """
<p class="exclusion-footnote">
*Reads mapped to an exclusion index (e.g. if mapping to an exclusion index containing
ribosomal RNA, this percentage is ribosomal RNA content).
</p>
"""

    if final_num_reads>0:
        rows.append(f"<tr><td>Final number of reads:</td><td>{final_num_reads:,}</td></tr>")

    table_html= "<table class='read-processing-table' style='width:60%; border-collapse:collapse; margin:0.8em 0;'>" + "\n".join(rows) + "</table>"

    html=[]
    html.append("<div class='read-processing-box' style='border:1px solid #888; background-color:#fafafa; padding:1em; margin:1em 0; border-radius:5px;'>")
    html.append("<h3>Read processing:</h3>")
    html.append(table_html)
    if ddd_note:
        html.append(ddd_note)
    if footnote_html:
        html.append(footnote_html)
    html.append("</div>")
    return "\n".join(html)


def generate_per_sample_alignment_html(sample_name, metrics, summary):
    """
    Reverting to "Number of input reads" (not Final number of reads).
    """
    if not metrics["input_reads"]:
        return ""
    inp= metrics.get("input_reads",0)
    pct_u= metrics.get("unique_pct",0.0)
    html=[]
    html.append("<div class='alignment-box' style='border:1px solid #339; background-color:#eef; padding:1em; margin:1em 0; border-radius:5px;'>")
    html.append("<h3>Alignment results (STARsolo):</h3>")

    html.append(f"""
    <table class="metric-table" style="width:50%; border-collapse:collapse; margin:0.8em 0;">
      <tr>
        <td style='border:1px solid #ccc; padding:8px 12px;'>Number of input reads</td>
        <td style='border:1px solid #ccc; padding:8px 12px;'><strong>{inp:,}</strong></td>
      </tr>
      <tr>
        <td style='border:1px solid #ccc; padding:8px 12px;'>Uniquely mapped</td>
        <td style='border:1px solid #ccc; padding:8px 12px;'><strong>{pct_u:.2f}%</strong></td>
      </tr>
      <tr>
        <td style='border:1px solid #ccc; padding:8px 12px;'>Uniquely mapped within genes (including introns)</td>
        <td style='border:1px solid #ccc; padding:8px 12px;'>{summary['map_genef']*100:.2f}%</td>
      </tr>
      <tr>
        <td style='border:1px solid #ccc; padding:8px 12px;'>Uniquely mapped within genes (excluding introns)</td>
        <td style='border:1px solid #ccc; padding:8px 12px;'>{summary['map_gene']*100:.2f}%</td>
      </tr>
      <tr>
        <td style='border:1px solid #ccc; padding:8px 12px;'>Sequencing saturation (1 - N<sub>umi</sub> / N<sub>reads</sub>)</td>
        <td style='border:1px solid #ccc; padding:8px 12px;'>{summary['sat']*100:.2f}%</td>
      </tr>
    </table>
    """)

    mm_mult= bullet_pct("Multimapped (multiple loci)", metrics.get("multimapped_pct",0.0), inp)
    mm_too= bullet_pct("Multimapped (too many loci)", metrics.get("multimapped_too_many_loci_pct",0.0), inp)
    un_mm= bullet_pct("Unmapped (too many mismatches)", metrics.get("unmapped_too_many_mismatches_pct",0.0), inp)
    un_sh= bullet_pct("Unmapped (too short)", metrics.get("unmapped_too_short_pct",0.0), inp)
    un_ot= bullet_pct("Unmapped (other)", metrics.get("unmapped_other_pct",0.0), inp, raw_count=metrics.get("unmapped_other_count",0))
    chm= bullet_pct("Chimeric reads", metrics.get("chimeric_pct",0.0), inp, raw_count=metrics.get("chimeric_count",0))

    html.append(f"""
<details>
<summary>Unmapped/Multimapped reads</summary>
<ul>
  <li>{mm_mult}</li>
  <li>{mm_too}</li>
  <li>{un_mm}</li>
  <li>{un_sh}</li>
  <li>{un_ot}</li>
  <li>{chm}</li>
</ul>
</details>
""")

    if sample_name!="merged":
        a_in= metrics.get("avg_input_read_length",0)
        a_map= metrics.get("avg_mapped_length",0)
        ns_tot= metrics.get("num_splices_total",0)
        ns_anno= metrics.get("num_splices_annotated",0)
        html.append(f"""
<details>
<summary>Other details</summary>
<ul>
<li>Average input read length: {a_in}</li>
<li>Average mapped length: {a_map}</li>
<li>Number of splices (total): {ns_tot:,}</li>
<li>Number of splices (annotated): {ns_anno:,}</li>
</ul>
</details>
""")

    html.append("</div>")
    return "\n".join(html)

def generate_merged_summary_table(per_sample_metrics, per_sample_fastqc, per_sample_kallisto):
    """
    Build a merged summary table with dynamic columns:
      - “Total reads”: STARsolo input_reads if >0, else kallisto n_processed if >0
      - “Uniquely mapped (%) – STARsolo”: only if any sample has unique_pct
      - “Pseudoalignment (%) – kallisto”: only if any sample has p_pseudoaligned

    Samples with neither STARsolo nor kallisto data for “Total reads” get a blank cell.
    """
    # Determine which columns to include
    include_star   = any(
        m.get("input_reads", 0) > 0 or m.get("unique_pct", None) is not None
        for m in per_sample_metrics.values() if m
    )
    if include_star:
        if per_sample_metrics[next(iter(per_sample_metrics))]["input_reads"] == 0:
            include_star = False
    include_kall   = any(
        k and k.get("p_pseudoaligned", None) is not None
        for k in per_sample_kallisto.values()
    )
    # If no STARsolo and no kallisto, nothing to show
    if not include_star and not include_kall:
        return "<p>No STARsolo or Kallisto data available for merged summary.</p>"

    # Build header
    cols = ["Sample"]
    if include_star or include_kall:
        cols.append("Total reads")
    if include_star:
        cols.append("Uniquely mapped (%) – STARsolo")
    if include_kall:
        cols.append("Pseudoalignment (%) – kallisto")

    # Start table
    html = []
    html.append("<table class='summary-table' style='border-collapse:collapse; margin-top:0.5em;'>")
    # header row
    html.append("<tr>" + "".join(
        f"<th style='border:1px solid #ccc; padding:8px;'>{col}</th>" for col in cols
    ) + "</tr>")

    # data rows (skip the 'merged' key itself)
    for sample, metrics in per_sample_metrics.items():
        if sample == "merged":
            continue
        row = []
        # Sample name
        row.append(f"<td style='border:1px solid #ccc; padding:8px;'>{sample}</td>")

        # Total reads: STARsolo input_reads else kallisto n_processed
        if include_star or include_kall:
            val = ""
            if metrics and metrics.get("input_reads", 0) > 0:
                val = f"{metrics['input_reads']:,}"
            else:
                kall = per_sample_kallisto.get(sample)
                if kall and kall.get("n_processed", 0) > 0:
                    val = f"{kall['n_processed']:,}"
            row.append(f"<td style='border:1px solid #ccc; padding:8px;'>{val}</td>")

        # Uniquely mapped (%) – STARsolo
        if include_star:
            pct = ""
            if metrics and metrics.get("unique_pct", None) is not None:
                pct = f"{metrics['unique_pct']:.2f}%"
            row.append(f"<td style='border:1px solid #ccc; padding:8px;'>{pct}</td>")

        # Pseudoalignment (%) – kallisto
        if include_kall:
            kpct = ""
            kall = per_sample_kallisto.get(sample)
            if kall and kall.get("p_pseudoaligned", None) is not None:
                kpct = f"{kall['p_pseudoaligned']:.1f}%"
            row.append(f"<td style='border:1px solid #ccc; padding:8px;'>{kpct}</td>")

        html.append("<tr>" + "".join(row) + "</tr>")

    html.append("</table>")
    return "\n".join(html)

###################################################################
def main(base_dir, run_name):
    pattern= os.path.join(base_dir,"results",f"output_*_starsolo_{run_name}")
    dirs= glob.glob(pattern)
    rx= re.compile(r"output_(.+)_starsolo_")
    all_samples=set()
    for d in dirs:
        mm= rx.search(os.path.basename(d))
        if mm:
            all_samples.add(mm.group(1))
    if len(all_samples) == 0:
        pattern= os.path.join(base_dir,"results",f"output_*_kallisto_{run_name}")
        dirs= glob.glob(pattern)
        rx= re.compile(r"output_(.+)_kallisto_")
        all_samples=set()
        for d in dirs:
            mm= rx.search(os.path.basename(d))
            if mm:
                all_samples.add(mm.group(1))
    all_samples= sorted(all_samples)

    per_sample_metrics={}
    per_sample_fastqc={}
    per_sample_exclusion={}
    per_sample_ligation={}
    per_sample_trimmed={}
    per_sample_cell_reads={}
    per_sample_anndata_geneFull={}
    per_sample_anndata_gene={}
    per_sample_anndata_kallisto_geneFull={}
    per_sample_anndata_kallisto_gene={}
    per_sample_kallisto = {}
    per_sample_pseudobulk = {}
    per_sample_summary = {}

    for sample_name in all_samples:
        star_log= os.path.join(base_dir,"results",f"output_{sample_name}_starsolo_{run_name}","Log.final.out")
        mm= parse_log_file(star_log)
        per_sample_metrics[sample_name]= mm

        fq_r1= os.path.join(base_dir,"qc",f"{sample_name}_R1_fastqc.html")
        # parse fastqc total
        fq_val= parse_fastqc_total_reads(fq_r1)
        per_sample_fastqc[sample_name]= fq_val

        ex_file= os.path.join(base_dir,f"{sample_name}.exclusion_stats.txt")
        ex_data= parse_exclusion_file(ex_file)
        per_sample_exclusion[sample_name]= ex_data

        lig_file= os.path.join(base_dir,f"{sample_name}.ligation_efficiency.txt")
        lig_data= parse_ligation_file(lig_file)
        per_sample_ligation[sample_name]= lig_data

        c_val= parse_trimmed_qc_files(base_dir, sample_name)
        per_sample_trimmed[sample_name]= c_val

        cell_reads_list = parse_cell_reads_stats(base_dir, sample_name, run_name)
        per_sample_cell_reads[sample_name] = cell_reads_list

        gf = parse_summary_csv(base_dir, sample_name, run_name, "GeneFull")
        g  = parse_summary_csv(base_dir, sample_name, run_name, "Gene")
        sat         = gf.get("Sequencing Saturation", g.get("Sequencing Saturation", 0.0))
        map_genef   = gf.get("Reads Mapped to GeneFull: Unique GeneFull", 0.0)
        map_gene    = g .get("Reads Mapped to Gene: Unique Gene",     0.0)
        per_sample_summary[sample_name] = {"sat": sat, "map_genef": map_genef, "map_gene": map_gene}



        per_sample_kallisto[sample_name] = parse_kallisto_info(base_dir, sample_name, run_name)

        gf_adata, g_adata = read_anndata_objects(base_dir, sample_name, run_name)
        per_sample_anndata_geneFull[sample_name] = gf_adata
        per_sample_anndata_gene[sample_name]     = g_adata
        gfk_adata, gk_adata = load_kallisto_anndata(base_dir, sample_name, run_name)
        per_sample_anndata_kallisto_geneFull[sample_name] = gfk_adata
        per_sample_anndata_kallisto_gene[sample_name] = gk_adata

        mtx_dir = os.path.join(base_dir, "results", f"output_{sample_name}_kallisto_{run_name}", "counts_unfiltered")
        cats = ["Nascent","Mature","Ambiguous"]
        paths = {c: os.path.join(mtx_dir, f"cells_x_genes.{c.lower()}.mtx") for c in cats}
        pseudo = {}
        pseudo_exists = False
        for c in cats:
            p = paths[c]
            if os.path.exists(p):
                tot, per_cell = parse_mtx_pseudobulk(p)
                pseudo_exists = True
            else:
                tot, per_cell = 0.0, {}
            pseudo[c] = (tot, per_cell)
        if pseudo_exists:
            per_sample_pseudobulk[sample_name] = pseudo
        else:
            per_sample_pseudobulk[sample_name] = None


    merged_geneFull = create_merged_anndata(per_sample_anndata_geneFull, all_samples)
    merged_gene = create_merged_anndata(per_sample_anndata_gene, all_samples)
    per_sample_anndata_geneFull["merged"] = merged_geneFull
    per_sample_anndata_gene["merged"] = merged_gene
    merged_kallisto_geneFull = create_merged_anndata(per_sample_anndata_kallisto_geneFull, all_samples)
    merged_kallisto_gene = create_merged_anndata(per_sample_anndata_kallisto_gene, all_samples)
    per_sample_anndata_kallisto_geneFull["merged"] = merged_kallisto_geneFull
    per_sample_anndata_kallisto_gene["merged"] = merged_kallisto_gene

    # Build merged Kallisto summary
    all_vals = [v for v in per_sample_kallisto.values() if v is not None]
    if all_vals:
        total_processed     = sum(v["n_processed"]     for v in all_vals)
        total_pseudoaligned = sum(v["n_pseudoaligned"] for v in all_vals)
        total_unique        = sum(v["n_unique"]        for v in all_vals)
        p_all   = (total_pseudoaligned / total_processed) * 100 if total_processed else 0.0
        p_uniq  = (total_unique        / total_processed) * 100 if total_processed else 0.0
        per_sample_kallisto["merged"] = {
            "n_processed":     total_processed,
            "p_pseudoaligned": p_all,
            "p_unique":        p_uniq,
            "call":            None
        }
    else:
        per_sample_kallisto["merged"] = None


    total_reads = sum(per_sample_metrics[s]["input_reads"] for s in all_samples)
    if total_reads != 0:
        sat_agg       = sum(per_sample_summary[s]["sat"]      * per_sample_metrics[s]["input_reads"] for s in all_samples) / total_reads
        map_genef_agg = sum(per_sample_summary[s]["map_genef"] * per_sample_metrics[s]["input_reads"] for s in all_samples) / total_reads
        map_gene_agg  = sum(per_sample_summary[s]["map_gene"]  * per_sample_metrics[s]["input_reads"] for s in all_samples) / total_reads
    else:
        sat_agg = None
        map_genef_agg = None
        map_gene_agg = None
    per_sample_summary["merged"] = {"sat": sat_agg, "map_genef": map_genef_agg, "map_gene": map_gene_agg}


    # build merged alignment
    merged_info={
        "input_reads": 0,
        "multimapped_count": 0,
        "multimapped_too_many_loci_count": 0,
        "unmapped_too_many_mismatches_count": 0,
        "unmapped_too_short_count": 0,
        "unmapped_other_count": 0,
        "chimeric_count": 0,
        "unique_count": 0,
        "avg_input_read_length_sum": 0.0,
        "avg_mapped_length_sum": 0.0,
        "num_splices_total_sum": 0,
        "num_splices_annotated_sum": 0,
        "sample_count": 0
    }
    for s,m in per_sample_metrics.items():
        n= m["input_reads"]
        merged_info["input_reads"]+= n
        ucount= int(round((m.get("unique_pct",0)/100.0)*n))
        merged_info["unique_count"]+= ucount

        mm_ct= int(round((m.get("multimapped_pct",0)/100.0)*n))
        merged_info["multimapped_count"]+= mm_ct

        mm_too= int(round((m.get("multimapped_too_many_loci_pct",0)/100.0)*n))
        merged_info["multimapped_too_many_loci_count"]+= mm_too

        unm_miss= int(round((m.get("unmapped_too_many_mismatches_pct",0)/100.0)*n))
        merged_info["unmapped_too_many_mismatches_count"]+= unm_miss

        unm_sh= int(round((m.get("unmapped_too_short_pct",0)/100.0)*n))
        merged_info["unmapped_too_short_count"]+= unm_sh

        merged_info["unmapped_other_count"]+= m["unmapped_other_count"]
        merged_info["chimeric_count"]+= m["chimeric_count"]
        merged_info["avg_input_read_length_sum"]+= m["avg_input_read_length"]
        merged_info["avg_mapped_length_sum"]+= m["avg_mapped_length"]
        merged_info["num_splices_total_sum"]+= m["num_splices_total"]
        merged_info["num_splices_annotated_sum"]+= m["num_splices_annotated"]
        merged_info["sample_count"]+=1

    sc= merged_info["sample_count"]
    tot_reads= merged_info["input_reads"]
    merged_metrics={}
    if tot_reads>0 and sc>0:
        merged_metrics["input_reads"]= tot_reads
        merged_metrics["unique_pct"]= (merged_info["unique_count"]/tot_reads)*100.0
        merged_metrics["multimapped_pct"]= (merged_info["multimapped_count"]/tot_reads)*100.0
        merged_metrics["multimapped_too_many_loci_pct"]= (merged_info["multimapped_too_many_loci_count"]/tot_reads)*100.0
        merged_metrics["unmapped_too_many_mismatches_pct"]= (merged_info["unmapped_too_many_mismatches_count"]/tot_reads)*100.0
        merged_metrics["unmapped_too_short_pct"]= (merged_info["unmapped_too_short_count"]/tot_reads)*100.0
        merged_metrics["unmapped_other_pct"]= (merged_info["unmapped_other_count"]/tot_reads)*100.0
        merged_metrics["unmapped_other_count"]= merged_info["unmapped_other_count"]
        merged_metrics["chimeric_pct"]= (merged_info["chimeric_count"]/tot_reads)*100.0
        merged_metrics["chimeric_count"]= merged_info["chimeric_count"]
        merged_metrics["avg_input_read_length"]= merged_info["avg_input_read_length_sum"]/sc
        merged_metrics["avg_mapped_length"]= merged_info["avg_mapped_length_sum"]/sc
        merged_metrics["num_splices_total"]= merged_info["num_splices_total_sum"]
        merged_metrics["num_splices_annotated"]= merged_info["num_splices_annotated_sum"]
    else:
        for k in ["input_reads","unique_pct","multimapped_pct","multimapped_too_many_loci_pct",
                  "unmapped_too_many_mismatches_pct","unmapped_too_short_pct","unmapped_other_pct","chimeric_pct"]:
            merged_metrics[k]=0
    per_sample_metrics["merged"]= merged_metrics

    all_cell_reads = []
    for s in all_samples:
        all_cell_reads.extend(per_sample_cell_reads[s])

    per_sample_cell_reads["merged"] = all_cell_reads

    # merge ligation
    def merge_ligs(L):
        mbc={}
        s_a=0
        for ld in L:
            for bc_c,vv in ld["barcode_dict"].items():
                mbc[bc_c]= mbc.get(bc_c,0)+ vv
            s_a+= ld["total_reads_post_trimming"]
        return {"barcode_dict":mbc, "total_reads_post_trimming": s_a}
    lig_list=[ per_sample_ligation[x] for x in all_samples]
    merged_lig= merge_ligs(lig_list)
    per_sample_ligation["merged"]= merged_lig

    # merge exclusion
    def merge_excs(E):
        te=0; tr=0
        for e in E:
            te+= e["excluded_count"]
            tr+= e["excluded_total"]
        mm={"excluded_count":te,"excluded_total":tr,"excluded_pct":0.0}
        if tr>0:
            mm["excluded_pct"]= (te/tr)*100.0
        return mm
    exc_list=[ per_sample_exclusion[x] for x in all_samples]
    merged_exc= merge_excs(exc_list)
    per_sample_exclusion["merged"]= merged_exc

    # merge fastqc
    tot_f= sum(per_sample_fastqc[x] for x in all_samples)
    per_sample_fastqc["merged"]= tot_f

    # merge trimmed
    tot_t= sum(per_sample_trimmed[x] for x in all_samples)
    per_sample_trimmed["merged"]= tot_t

    ############################################
    # Build the tabs => "merged" + each sample
    ############################################
    tab_order= ["merged"] + list(all_samples)

    tab_links=[]
    tab_contents=[]

    for i,sample in enumerate(tab_order):
        tab_id= f"tab_{sample}"
        class_act= "active" if i==0 else ""
        tab_links.append(f'<button class="tablinks {class_act}" onclick="openTab(event, \'{tab_id}\')">{sample}</button>')

        # read-processing
        c_val= per_sample_trimmed.get(sample,0)
        lig_dat= per_sample_ligation.get(sample,{})
        exc_dat= per_sample_exclusion.get(sample,{})
        fq_val= per_sample_fastqc.get(sample,0)
        final_reads= per_sample_metrics.get(sample,{}).get("input_reads",0)

        # alignment
        alg_dat= per_sample_metrics.get(sample,{})
        align_html= generate_per_sample_alignment_html(sample, alg_dat, per_sample_summary[sample])
        kall = per_sample_kallisto.get(sample)
        if sample != "merged":
            div_pseudomtx = create_nascent_mature_ambiguous_bars(per_sample_pseudobulk[sample], sample)
        else:
            merged = { c:(0, {}) for c in ["Nascent","Mature","Ambiguous","Total"] }
            pseudo_exists = False
            for s in all_samples:
                if not per_sample_pseudobulk[s]:
                    continue
                else:
                    pseudo_exists = True
                for c,(tot,pc) in per_sample_pseudobulk[s].items():
                    m_tot, m_pc = merged[c]
                    m_tot += tot
                    for cell,val in pc.items():
                        m_pc[cell] = m_pc.get(cell,0) + val
                    merged[c] = (m_tot, m_pc)
            if not pseudo_exists:
                merged = None
            div_pseudomtx = create_nascent_mature_ambiguous_bars(merged, "merged")
        if kall:
            kallisto_html = f"""
        <div class="alignment-box" style="border:1px solid #339; background-color:#eef; padding:1em; margin:1em 0; border-radius:5px;">
          <h3>Kallisto pseudoalignment:</h3>
          <table class="metric-table" style="width:50%; border-collapse:collapse; margin:0.8em 0;">
            <tr>
              <td style="border:1px solid #ccc; padding:8px 12px;">Number of input reads</td>
              <td style="border:1px solid #ccc; padding:8px 12px;"><strong>{kall['n_processed']:,}</strong></td>
            </tr>
            <tr>
              <td style="border:1px solid #ccc; padding:8px 12px;">Percent pseudoaligned</td>
              <td style="border:1px solid #ccc; padding:8px 12px;"><strong>{kall['p_pseudoaligned']:.1f}%</strong></td>
            </tr>
            <tr>
              <td style="border:1px solid #ccc; padding:8px 12px;">Percent pseudoaligned to unique transcripts</td>
              <td style="border:1px solid #ccc; padding:8px 12px;"><strong>{kall['p_unique']:.1f}%</strong></td>
            </tr>
          </table>
            """
            if kall['call']:
                kallisto_html += (f"""
                 <details>
                   <summary style="cursor: pointer;">Show kallisto run details</summary>
                   <pre style="white-space: pre-wrap; word-wrap: break-word; background:#f7f7f7; padding:0.5em; border-radius:4px; margin-top:0.5em;"><code>{kall['call']}</code></pre>
                 </details>
               """)
            kallisto_html += div_pseudomtx
            kallisto_html += ("</div>")
        else:
            kallisto_html = ""

        if not final_reads and kall and kall['n_processed']:
            final_reads = kall['n_processed']

        # read processing (add to the HTML)
        rp_html= generate_read_processing_html(fq_val,c_val,lig_dat,exc_dat,final_reads)

        bc_html = "" # Barcode distribution
        barcode_box = create_suffix_umi_comparison_box(per_sample_anndata_geneFull.get(sample),
                "STARsolo",
                per_sample_anndata_kallisto_geneFull.get(sample),
                "kallisto",
                sample_name=sample
        )
        if barcode_box:
            bc_html += "<div class='read-processing-box' style='border:1px solid #888;"
            bc_html += "background-color:#fafafa;padding:1em;margin:1em 0;border-radius:5px;'>"
            bc_html += "<h3>Barcode representation</h3>"
            bc_html += "<div style='max-height:360px; max-width: 1000px; background-color: #ffffff; overflow-y:auto; padding-right:10px; border: 1px solid #888'>"
            bc_html += barcode_box
            bc_html += "</div>"
            bc_html += "</div>"


        # fastqc toggles
        if sample=="merged":
            fc_html= ""
        else:
            fc_html= generate_fastqc_iframe_toggles(sample, base_dir)

        sat_div = create_sequencing_saturation_plot(
            per_sample_cell_reads[sample],
            sample_name=sample,
            is_merged=(sample=="merged")
        )

        scatter_html = create_mouse_human_scatter(
            adata_geneFull=per_sample_anndata_geneFull[sample] if per_sample_anndata_geneFull[sample] else per_sample_anndata_kallisto_geneFull[sample],
            sample_name=sample,
            label_="STARsolo" if per_sample_anndata_geneFull[sample] else "kallisto",
            is_merged=(sample=="merged")
        )

        knee_div = create_knee_plot_html_from_anndata(
            per_sample_anndata_geneFull[sample],
            per_sample_anndata_gene[sample],
            per_sample_anndata_kallisto_geneFull[sample],
            per_sample_anndata_kallisto_gene[sample],
            sample_name=sample,
            is_merged=(sample=="merged")
        )


        # analysis => knee plot
        if sample=="merged":
            # pass all samples => unify columns
            # knee_div= create_knee_plot_html(all_samples, base_dir, run_name, is_merged=True)
            analysis_box=[]
            analysis_box.append("<div class='analysis-box' style='border:1px solid #888; background-color:#ffffff; padding:1em; margin:1em 0; border-radius:5px;'>")
            analysis_box.append("<h3>Analysis</h3>")
            analysis_box.append(knee_div)
            analysis_box.append(sat_div)
            analysis_box.append(scatter_html)
            analysis_box.append("</div>")

            # separate box => per-sample summary
            sum_table= generate_merged_summary_table(per_sample_metrics, per_sample_fastqc, per_sample_kallisto)
            summary_box=[]
            summary_box.append("<div class='analysis-box' style='border:1px solid #888; background-color:#fafafa; padding:1em; margin:1em 0; border-radius:5px;'>")
            summary_box.append("<h3>Per-sample summary</h3>")
            summary_box.append(sum_table)
            summary_box.append("</div>")

            extra_html= "\n".join(analysis_box) + "\n" + "\n".join(summary_box)
        else:
            #knee_div= create_knee_plot_html([sample], base_dir, run_name, is_merged=False)
            analysis_box=[]
            analysis_box.append("<div class='analysis-box' style='border:1px solid #888; background-color:#ffffff; padding:1em; margin:1em 0; border-radius:5px;'>")
            analysis_box.append("<h3>Analysis</h3>")
            analysis_box.append(knee_div)
            analysis_box.append(sat_div)
            analysis_box.append(scatter_html)
            analysis_box.append("</div>")
            extra_html= "\n".join(analysis_box)

        content_html= rp_html + "\n" + align_html + "\n" + kallisto_html + "\n" + fc_html + "\n" + extra_html + "\n" + bc_html

        style_disp= "block" if i==0 else "none"
        tab_contents.append(f'<div id="{tab_id}" class="tabcontent" style="display:{style_disp};">\n{content_html}\n</div>')

    ############################################
    # Final HTML
    ############################################
    final_html=[]
    final_html.append("<html><head><meta charset='UTF-8'><title>SWIFT-seq report</title>")
    final_html.append(f"<script type='text/javascript'>{_PLOTLY_JS}</script>")
    final_html.append("""
<style>
body {
  font-family: Arial, sans-serif;
  margin: 20px;
}

details summary {
  cursor: pointer;
}

.tab {
  border-bottom: 1px solid #ccc;
  margin-bottom: 0.5em;
}
.tab button {
  background-color: #f1f1f1;
  border: none;
  outline: none;
  cursor: pointer;
  padding: 8px 16px;
  margin-right: 4px;
}
.tab button.active {
  background-color: #ccc;
}
.tabcontent {
  display: none;
  padding: 0.5em 0;
}

/* old table style */
.read-processing-box {
  border:1px solid #888;
  background-color:#fafafa;
  padding:1em;
  margin:1em 0;
  border-radius:5px;
}
.read-processing-table {
  border-collapse: collapse;
  margin: 0.8em 0;
  width: 60%;
}
.read-processing-table td {
  border: 1px solid #ccc;
  padding: 8px 12px;
  vertical-align: middle;
}

.alignment-box {
  border: 1px solid #339;
  background-color: #eef;
  padding: 1em;
  margin: 1em 0;
  border-radius: 5px;
}
.metric-table {
  border-collapse: collapse;
  margin: 0.8em 0;
  width: 50%;
}
.metric-table td {
  border: 1px solid #ccc;
  padding: 8px 12px;
  vertical-align: middle;
}

.summary-table {
  border-collapse: collapse;
  margin-top: 0.5em;
}
.summary-table th, .summary-table td {
  border: 1px solid #ccc;
  padding: 8px;
}

.analysis-box {
  border:1px solid #888;
  background-color:#fafafa;
  padding:1em;
  margin:1em 0;
  border-radius:5px;
}

.fastqc-toggles {
  margin:1em 0;
  padding:0;
}
</style>

<script>
function openTab(evt, tabId){
  var tabcontent = document.getElementsByClassName("tabcontent");
  for(var i=0; i<tabcontent.length; i++){
    tabcontent[i].style.display= "none";
  }
  var tablinks= document.getElementsByClassName("tablinks");
  for(var i=0; i<tablinks.length; i++){
    tablinks[i].className= tablinks[i].className.replace(" active","");
  }
  document.getElementById(tabId).style.display= "block";
  evt.currentTarget.className+= " active";
}

// toggling for the R1/R2 iFrame
function toggleFastQC(divID){
  var el = document.getElementById(divID);
  if(!el) return;
  if(el.style.display==="none"){
    el.style.display="block";
  } else {
    el.style.display="none";
  }
}
</script>
""")
    final_html.append("</head><body>")
    final_html.append("<h1>SWIFT-seq report</h1>")
    final_html.append('<div class="tab">')
    final_html.append("".join(tab_links))
    final_html.append("</div>")
    final_html.append("".join(tab_contents))
    final_html.append("</body></html>")

    out_path= os.path.join(base_dir,"results",f"report_{run_name}.html")
    with open(out_path,"w",encoding="utf-8") as f:
        f.write("\n".join(final_html))

    print(f"Report generated: {out_path}")


if __name__=="__main__":
    if len(sys.argv)<3:
        print(f"Usage: {sys.argv[0]} <base_dir> <run_name>")
        sys.exit(1)
    base_dir= sys.argv[1]
    run_name= sys.argv[2]
    main(base_dir, run_name)

