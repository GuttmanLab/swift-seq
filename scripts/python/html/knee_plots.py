# knee_plots.py
import os
import scipy.sparse as sp
from scipy.io import mmread
import numpy as np
import re
from collections import defaultdict

import plotly.graph_objs as go
import plotly.offline as pyo

_PLOTLY_JS_EMBEDDED = False          # global flag
_PLOTLY_JS = pyo.get_plotlyjs()          # full minified plotly.js
_PLOTLY_READY = False                    # flag toggled after we inject JS

def plot_div_once(fig):
    """
    Return an HTML <div> for a Plotly figure.
    The very first time this is called it embeds full plotly.js
    (include_plotlyjs=True).  Every subsequent call sets
    include_plotlyjs=False so the JS library isn’t duplicated.
    """
    global _PLOTLY_JS_EMBEDDED

    div = pyo.plot(
        fig,
        include_plotlyjs=False, # not _PLOTLY_JS_EMBEDDED,   # True only the 1st time
        output_type="div",
        show_link=False
    )
    _PLOTLY_JS_EMBEDDED = True      # flip the switch after first call
    return div

def _compute_knee_xy(umis_per_cell):
    """
    For each unique UMI value u, figure out how many barcodes have >= u.
    We'll produce (x, y) for a knee plot, x ascending, y descending.
    """
    from collections import Counter
    sums_int = umis_per_cell.astype(int)
    c = Counter(sums_int)
    items = sorted(c.items())  # ascending by UMIs
    total_cells = len(sums_int)

    # partial sums
    import itertools
    from itertools import accumulate
    counts_only = [it[1] for it in items]
    sums_cumu = list(accumulate(counts_only))  # ascending cumsum

    x_vals = []
    y_vals = []
    for i, (umi_val, cnt_here) in enumerate(items):
        if i == 0:
            less_count = 0
        else:
            less_count = sums_cumu[i - 1]
        at_least = total_cells - less_count
        x_vals.append(umi_val)
        y_vals.append(at_least)

    return x_vals, y_vals

def _load_concat_matrices(sample_list, base_dir, run_name, which="GeneFull"):
    """
    Memory-efficient horizontal stack (columns) of each sample's matrix:
      base_dir/results/output_<sample>_starsolo_<run_name>/Solo.out/<which>/raw/matrix.mtx
    Return a CSR or None if not found for that 'which' in any sample.
    """
    sp_mats = []
    gene_dim = None

    for s in sample_list:
        mat_path = os.path.join(
            base_dir, "results", f"output_{s}_starsolo_{run_name}",
            "Solo.out", which, "raw", "matrix.mtx"
        )
        if not os.path.exists(mat_path):
            continue
        M_coo = mmread(mat_path).tocsr()  # still sparse
        if gene_dim is None:
            gene_dim = M_coo.shape[0]
        else:
            if M_coo.shape[0] != gene_dim:
                print(f"WARNING: mismatch in gene dimension for sample {s}.")
        sp_mats.append(M_coo)

    if len(sp_mats) == 0:
        return None

    big_mat = sp.hstack(sp_mats, format="csr")
    return big_mat

def create_knee_plot_html_from_anndata(
    adata_geneFull,
    adata_gene,
    adata_kallisto_geneFull,
    adata_kallisto_gene,
    sample_name,
    is_merged=False
):
    """
    Creates a two-line knee plot (GeneFull vs. Gene) using the logic:
      knee = np.sort(adata.X.sum(axis=1))[::-1]
      x = knee
      y = range(len(knee))

    - 'Genes including introns' is visible by default
    - 'Genes excluding introns' is hidden by default
    - Uses log–log axes, line width=5, bigger fonts, axis lines, 1000×500 figure

    If adata_geneFull or adata_gene is None, we skip that line.
    Returns an offline Plotly <div> with the figure.
    """

    data_traces = []

    # Genes including introns
    if adata_geneFull:
        x_incl, y_incl = _compute_knee_xy(np.array(adata_geneFull.X.sum(axis=1)).ravel())
        text_incl = [f"UMIs: {x}\nBarcodes: {y}" for x, y in zip(x_incl, y_incl)]
        trace_incl = go.Scatter(
            x=x_incl,
            y=y_incl,
            mode="lines+markers",
            name="Genes including introns (STARsolo)",
            text=text_incl,
            hoverinfo="text",
            line=dict(width=3),
            visible=True  # shown by default
        )
        data_traces.append(trace_incl)

    # Genes excluding introns
    if adata_gene:
        x_excl, y_excl = _compute_knee_xy(np.array(adata_gene.X.sum(axis=1)).ravel())
        text_excl = [f"UMIs: {x}\nBarcodes: {y}" for x, y in zip(x_excl, y_excl)]
        trace_excl = go.Scatter(
            x=x_excl,
            y=y_excl,
            mode="lines+markers",
            name="Genes excluding introns (STARsolo)",
            text=text_excl,
            hoverinfo="text",
            line=dict(width=3),
            visible="legendonly"  # hidden by default
        )
        data_traces.append(trace_excl)

    if adata_kallisto_geneFull:
        x_incl_k, y_incl_k = _compute_knee_xy(np.array(adata_kallisto_geneFull.X.sum(axis=1)).ravel())
        text_incl_k = [f"UMIs: {x}\nBarcodes: {y}" for x, y in zip(x_incl_k, y_incl_k)]
        trace_incl_k = go.Scatter(
            x=x_incl_k,
            y=y_incl_k,
            mode="lines+markers",
            name="Genes including introns (kallisto)",
            text=text_incl_k,
            hoverinfo="text",
            line=dict(width=3),
            visible="legendonly" if adata_geneFull else True
        )
        data_traces.append(trace_incl_k)
    if adata_kallisto_gene:
        x_excl_k, y_excl_k = _compute_knee_xy(np.array(adata_kallisto_gene.X.sum(axis=1)).ravel())
        text_excl_k = [f"UMIs: {x}\nBarcodes: {y}" for x, y in zip(x_excl_k, y_excl_k)]
        trace_excl_k = go.Scatter(
            x=x_excl_k,
            y=y_excl_k,
            mode="lines+markers",
            name="Genes excluding introns (kallisto)",
            text=text_excl_k,
            hoverinfo="text",
            line=dict(width=3),
            visible="legendonly"
        )
        data_traces.append(trace_excl_k)


    # If neither AnnData was provided, or both empty, return a note
    if not data_traces:
        name_str = "merged" if is_merged else sample_name
        return f"<p>No knee plot data for {name_str}.</p>"

    plot_title = "Knee plot"
    layout = go.Layout(
        # The "Analysis" text is not used here, since we want the box title instead
        title=plot_title,
        width=1000,
        height=500,
        showlegend=True,
        legend=dict(font=dict(size=16)),
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(
            type="log",
            title=dict(text="UMI counts",font=dict(size=18)),
            dtick=1,  # only powers of 10
            tickfont=dict(size=16),
            zeroline=False,
            showline=True,
            linecolor="black",
            linewidth=2
        ),
        yaxis=dict(
            type="log",
            title=dict(text="Number of barcodes",font=dict(size=18)),
            dtick=1,
            tickfont=dict(size=16),
            zeroline=False,
            showline=True,
            linecolor="black",
            linewidth=2
        ),
        hovermode="closest"
    )

    fig = go.Figure(data=data_traces, layout=layout)

    # produce offline HTML
    plot_div = plot_div_once(fig)
    #plot_div = pyo.plot(
    ##    fig,
    #    include_plotlyjs=True,  # offline mode => embed entire JS
    #    output_type="div",
    #    show_link=False
    #)
    return plot_div




def create_knee_plot_html(sample_names, base_dir, run_name, is_merged=False):
    """
    Build an offline Plotly figure (embedding entire plotly.js),
    that shows up to two lines:
      - Genes including introns (GeneFull) => visible
      - Genes excluding introns (Gene) => legendonly

    width=1000, height=500,
    log-log axes,
    x-axis, y-axis lines are visible (like a "frame" around the plot),
    bigger fonts (tick=16, title=18, legend=16).

    Return an HTML <div> with the figure + the needed JS for offline usage.
    """
    big_incl = _load_concat_matrices(sample_names, base_dir, run_name, which="GeneFull")
    big_excl = _load_concat_matrices(sample_names, base_dir, run_name, which="Gene")

    data_traces = []

    # Genes including introns
    if big_incl is not None:
        sums_incl = np.array(big_incl.sum(axis=0)).ravel()
        x_incl, y_incl = _compute_knee_xy(sums_incl)
        text_incl = [f"UMIs: {x}\nBarcodes: {y}" for x, y in zip(x_incl, y_incl)]
        trace_incl = go.Scatter(
            x=x_incl,
            y=y_incl,
            mode="lines+markers",
            name="Genes including introns (STARsolo)",
            text=text_incl,
            hoverinfo="text",
            line=dict(width=3),
            visible=True  # shown by default
        )
        data_traces.append(trace_incl)

    # Genes excluding introns
    if big_excl is not None:
        sums_excl = np.array(big_excl.sum(axis=0)).ravel()
        x_excl, y_excl = _compute_knee_xy(sums_excl)
        text_excl = [f"UMIs: {x}\nBarcodes: {y}" for x, y in zip(x_excl, y_excl)]
        trace_excl = go.Scatter(
            x=x_excl,
            y=y_excl,
            mode="lines+markers",
            name="Genes excluding introns (STARsolo)",
            text=text_excl,
            hoverinfo="text",
            line=dict(width=3),
            visible="legendonly"  # hidden by default
        )
        data_traces.append(trace_excl)

    plot_title = "Knee plot"
    layout = go.Layout(
        # The "Analysis" text is not used here, since we want the box title instead
        title=plot_title,
        width=1000,
        height=500,
        showlegend=True,
        legend=dict(font=dict(size=16)),
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(
            type="log",
            title=dict(text="UMI counts",font=dict(size=18)),
            dtick=1,  # only powers of 10
            tickfont=dict(size=16),
            zeroline=False,
            showline=True,
            linecolor="black",
            linewidth=2
        ),
        yaxis=dict(
            type="log",
            title=dict(text="Number of barcodes",font=dict(size=18)),
            dtick=1,
            tickfont=dict(size=16),
            zeroline=False,
            showline=True,
            linecolor="black",
            linewidth=2
        ),
        hovermode="closest"
    )

    fig = go.Figure(data=data_traces, layout=layout)

    # produce offline HTML
    plot_div = plot_div_once(fig)
    #plot_div = pyo.plot(
    #    fig,
    #    include_plotlyjs=True,  # offline mode => embed entire JS
    #    output_type="div",
    #    show_link=False
    #)
    return plot_div

def create_sequencing_saturation_plot(cell_reads_list, sample_name, is_merged=False):
    """
    Bins the cells by read count as follows:
       - 0..1000 in increments of 200
       - 1000..10000 in increments of 1000
       - 10000..100000 in increments of 5000
       - 100000..1000000 in increments of 25000
    For each bin, computes mean reads (X) and mean UMIs (Y).
    Then plots these bin means as a lines+markers scatter.

    cell_reads_list: list of (reads, umis) per cell
    sample_name:     name of this sample
    is_merged:       whether it's a merged sample (for the title)

    Returns an HTML <div> with a Plotly figure embedded (offline).
    """

    import plotly.graph_objs as go
    import plotly.offline as pyo

    # If no data
    if not cell_reads_list:
        name_str = "merged" if is_merged else sample_name
        return ""

    # We'll build bin edges manually:
    bin_edges = []

    # 0..1000 in increments of 250
    current = 0
    while current <= 1000:
        bin_edges.append(current)
        current += 250

    # 1000..10000 in increments of 1000
    edge = 1000
    while edge <= 10000:
        bin_edges.append(edge)
        edge += 1000

    # 10000..100000 in increments of 6000
    edge = 10000
    while edge <= 100000:
        bin_edges.append(edge)
        edge += 6000

    # 100000..1000000 in increments of 10000
    edge = 100000
    while edge <= 1000000:
        bin_edges.append(edge)
        edge += 10000

    # 1000000..10000000 in increments of 100000
    edge = 1000000
    while edge <= 10000000:
        bin_edges.append(edge)
        edge += 100000

    # 10000000..100000000 in increments of 1000000
    edge = 10000000
    while edge <= 100000000:
        bin_edges.append(edge)
        edge += 1000000

    # 100000000..1000000000 in increments of 10000000
    edge = 100000000
    while edge <= 1000000000:
        bin_edges.append(edge)
        edge += 10000000

    # 1000000000..10000000000 in increments of 100000000
    edge = 1000000000
    while edge <= 10000000000:
        bin_edges.append(edge)
        edge += 100000000

    # 10000000000..100000000000 in increments of 1000000000
    edge = 10000000000
    while edge <= 100000000000:
        bin_edges.append(edge)
        edge += 1000000000


    # Sort & unique bin_edges in case there's overlap
    bin_edges = sorted(set(bin_edges))

    # Sort cells by reads
    cells_sorted = sorted(cell_reads_list, key=lambda x: x[0])  # (reads, umis)

    x_means = []
    y_means = []

    # We'll iterate over consecutive bin pairs
    for i in range(len(bin_edges) - 1):
        bin_min = bin_edges[i]
        bin_max = bin_edges[i+1]

        # gather cells with reads in [bin_min, bin_max)
        # You could also do <= bin_max if you prefer inclusive, but let's do [min, max)
        chunk = [(r,u) for (r,u) in cells_sorted if (r >= bin_min and r < bin_max)]
        if not chunk:
            continue  # skip empty bin

        # compute average reads & UMIs in this bin
        sum_r = 0.0
        sum_u = 0.0
        for (reads_i, umis_i) in chunk:
            sum_r += reads_i
            sum_u += umis_i
        count = len(chunk)
        mean_reads = sum_r / count
        mean_umis  = sum_u / count

        x_means.append(mean_reads)
        y_means.append(mean_umis)

    # If for some reason we never got any bins, return a note
    if not x_means:
        name_str = "merged" if is_merged else sample_name
        return f"<p>No cells fall into any bins for {name_str}.</p>"

    # Build a lines+markers scatter
    scatter_trace = go.Scatter(
        x=x_means,
        y=y_means,
        mode="lines+markers",
        name="Mean Reads vs Mean UMIs (Custom binning)",
        hoverinfo="x+y",
        line=dict(width=3)
    )

    plot_title = f"Sequencing Saturation: {'merged' if is_merged else sample_name}"
    layout = go.Layout(
        title=plot_title,
        width=600,
        height=500,
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(
            title=dict(text="Mean reads per cell",font=dict(size=18)),
            tickfont=dict(size=16),
            showline=True,
            linecolor="black",
            linewidth=2
        ),
        yaxis=dict(
            title=dict(text="Mean unique molecules",font=dict(size=18)),
            tickfont=dict(size=16),
            showline=True,
            linecolor="black",
            linewidth=2
        ),
    )
    fig = go.Figure(data=[scatter_trace], layout=layout)


    #div_html = pyo.plot(fig, include_plotlyjs=True, output_type="div")
    div_html = plot_div_once(fig)
    return div_html

def create_mouse_human_scatter(
    adata_geneFull,
    sample_name,
    label_,
    is_merged=False
):
    """
    1) Check if AnnData has >= 50 genes starting 'ENSG' (human) and >= 50 genes
       starting 'ENSMUSG' (mouse) in adata.var['gene_id'] -> if not, return note.
    2) If yes => for each cell, sum UMIs from human genes, sum UMIs from mouse genes.
    3) X = mouse_UMI, Y = human_UMI (both are sums across relevant gene IDs).
    4) Classify each cell color: 'human' if ratio >=0.95, 'mouse' if ratio <=0.05, else 'gray'.
    5) Make scatter plot with x=mouse, y=human, color-coded points.

    Returns Plotly offline <div> for embedding in HTML.
    """
    import numpy as np
    import plotly.graph_objs as go
    import plotly.offline as pyo

    if adata_geneFull is None or adata_geneFull.n_obs == 0 or adata_geneFull.n_vars == 0:
        name_str = "merged" if is_merged else sample_name
        return "" # f"<p>No valid data for {name_str} in adata_geneFull.</p>"

    # We assume the gene IDs are in adata_geneFull.var['gene_id']
    if "gene_id" not in adata_geneFull.var.columns:
        return "" # f"<p>No 'gene_id' column in var. Can't check for human/mouse genes.</p>"

    gene_ids = adata_geneFull.var["gene_id"].values.astype(str)

    # Count how many start with 'ENSG' vs 'ENSMUSG'
    is_human_gene = np.char.startswith(gene_ids, "ENSG")
    is_mouse_gene = np.char.startswith(gene_ids, "ENSMUSG")

    n_human_genes = np.sum(is_human_gene)
    n_mouse_genes = np.sum(is_mouse_gene)

    # 1) Check the threshold => 50 each
    if n_human_genes < 50 or n_mouse_genes < 50:
        name_str = "merged" if is_merged else sample_name
        return ""
        # return f"<p>Not a human-mouse mixing experiment for {name_str} (need >=50 each, found {n_human_genes} human, {n_mouse_genes} mouse).</p>"

    # 2) For each cell, sum UMIs from human & mouse genes
    # We'll separate the columns for human vs. mouse genes:
    # If .X is huge, this might be memory-intensive, but it's the simplest direct approach
    # We'll do a partial approach: sum along axis=1 but only the relevant columns
    # watch out for sparse => do np.asarray(...) if needed

    # Indices for human vs. mouse genes
    human_cols = np.where(is_human_gene)[0]
    mouse_cols = np.where(is_mouse_gene)[0]

    # Summation:
    # if .X is e.g. csc, you can do .X[:, human_cols].sum(axis=1)
    # We'll do the simplest approach:
    human_sums = np.array(adata_geneFull.X[:, human_cols].sum(axis=1)).ravel()
    mouse_sums = np.array(adata_geneFull.X[:, mouse_cols].sum(axis=1)).ravel()
    tot_umis = human_sums + mouse_sums
    mask     = tot_umis >= 5 # boolean array
    human_sums = human_sums[mask]
    mouse_sums = mouse_sums[mask]

    # 3) For color classification:
    # ratio = human_sums / (human_sums + mouse_sums)
    # if ratio >=0.95 => 'human'
    # if ratio <=0.05 => 'mouse'
    # else => 'gray'
    total_sums = human_sums + mouse_sums
    eps = 1e-9
    ratio = np.divide(human_sums, total_sums + eps)  # prevent /0

    colors = []
    for r in ratio:
        if r >= 0.95:
            colors.append("red")   # or "human"
        elif r <= 0.05:
            colors.append("blue") # or "mouse"
        else:
            colors.append("gray")

    # We'll store them in a single trace => mode="markers"
    # x=mouse_sums, y=human_sums
    scatter_trace = go.Scattergl(
        x=mouse_sums + 0.1,
        y=human_sums + 0.1,
        mode="markers",
        marker=dict(
            size=6,
            color=colors
        ),
        name="Cells: Human vs Mouse mixing",
        hoverinfo="skip"
    )

    merged_str = "merged" if is_merged else sample_name
    layout = go.Layout(
        title=f"Human-Mouse mixing: {merged_str} ({label_})",
        width=800,
        height=600,
        xaxis=dict(
            title="Mouse UMI counts",
            zeroline=False,
            showline=True,
            linecolor="black",
            linewidth=2
        ),
        yaxis=dict(
            title="Human UMI counts",
            zeroline=False,
            showline=True,
            linecolor="black",
            linewidth=2
        )
    )

    updatemenus = [
    dict(
        type="buttons",
        direction="up",
        buttons=list([
            dict(
                args=[{"xaxis.type": "linear", "yaxis.type": "linear"}],
                label="Linear Scale",
                method="relayout"
            ),
            dict(
                args=[{"xaxis.type": "log", "yaxis.type": "log"}],
                label="Log Scale",
                method="relayout"
            )
        ])
    ),
    ]

    fig = go.Figure(data=[scatter_trace], layout=layout)
    fig.update_layout(updatemenus=updatemenus)
    # Return offline HTML <div>
    #scatter_div = pyo.plot(fig, include_plotlyjs=True, output_type="div")
    scatter_div = plot_div_once(fig)
    return scatter_div


def create_nascent_mature_ambiguous_bars(pseudobulk_dict, sample_name):
    """
    pseudobulk_dict maps each of "Nascent","Mature","Ambiguous" → (tot, per_cell_dict).
    We derive total_per_cell by summing those three per_cell dicts.
    """
    import numpy as np
    from collections import defaultdict

    if not pseudobulk_dict:
        return ""

    cats = ["Nascent","Mature","Ambiguous"]
    cat_labels = ["Nascent RNA", "Mature RNA", "Ambiguous (exon-only)"]

    # Extract totals and per-cell dicts
    totals   = {c: pseudobulk_dict[c][0] for c in cats}
    per_cell = {c: pseudobulk_dict[c][1] for c in cats}

    # Build total_per_cell by summing over the three categories
    total_per_cell = defaultdict(float)
    for c in cats:
        for cell, val in per_cell[c].items():
            total_per_cell[cell] += val

    # Find which cells have ≥500 UMIs total
    cells_500 = {cell for cell, tsum in total_per_cell.items() if tsum >= 500}

    # Compute filtered sums for each category
    filtered = {}
    for c in cats:
        s = 0.0
        for cell, val in per_cell[c].items():
            if cell in cells_500:
                s += val
        filtered[c] = s


    all_vals = [totals[c] for c in cats]
    filt_vals = [filtered[c] for c in cats]

    # compute percentages
    sum_all  = sum(all_vals)
    sum_filt = sum(filt_vals)

    all_pcts  = [f"{(v/sum_all*100):.1f}%" for v in all_vals]
    filt_pcts = [f"{(v/sum_filt*100):.1f}%" for v in filt_vals]


    # Build Plotly traces
    trace_all = go.Bar(
        x=[totals[c] for c in cats],
        y=cat_labels,
        orientation="h",
        name="All cells",
        marker=dict(color="steelblue"),
        text=all_pcts,
        textposition="inside",
        hovertemplate=("%{text}")
    )
    trace_filt = go.Bar(
        x=[filtered[c] for c in cats],
        y=cat_labels,
        orientation="h",
        name="Cells ≥500 UMIs",
        marker=dict(color="steelblue"),
        visible=False,
        text=filt_pcts,
        textposition="inside",
        hovertemplate=("%{text}")
    )


    # Toggle buttons
    updatemenu = [{
        "buttons": [
            {"method":"update","label":"All cells",
             "args":[{"visible":[True,False]},{"title":f"{sample_name}: All cells"}]},
            {"method":"update","label":"Cells ≥500 UMIs",
             "args":[{"visible":[False,True]},{"title":f"{sample_name}: Cells ≥500 UMIs"}]}
        ],
        "direction":"right","x":0,"y":1.2,"xanchor":"left","yanchor":"top"
    }]

    layout = go.Layout(
        barmode="overlay",
        width=700, height=315,
        updatemenus=updatemenu,
        xaxis=dict(title="Total counts"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=18)),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        font=dict(size=18)
    )

    fig = go.Figure(data=[trace_all, trace_filt], layout=layout)
    return plot_div_once(fig) # pyo.plot(fig, include_plotlyjs=False, output_type="div", show_link=False)


import re
import numpy as np
import plotly.graph_objs as go
import plotly.offline as pyo
from collections import defaultdict

def create_suffix_umi_comparison_box(
    adata1, label1,
    adata2=None, label2=None,
    sample_name=None
):
    """
    adata1 + label1: first AnnData + its name ("GeneFull")
    adata2 + label2: optional second AnnData + its label ("kallisto GeneFull")
    sample_name: used in the box title
    """
    # Helper to extract suffix→UMI sums
    def suffix_sums(adata):
        ids = adata.obs["id"].astype(str).values
        pat = re.compile(r"\[([^]]+)\]$")
        suf = [ (m.group(1) if (m:=pat.search(v)) else "__NA__")
                for v in ids ]
        umi = np.asarray(adata.X.sum(axis=1)).ravel()
        d = defaultdict(float)
        for s,u in zip(suf, umi):
            d[s] += u
        return d

    # Collect suffix sums for each dataset
    sums1 = suffix_sums(adata1) if adata1 is not None else {}
    sums2 = suffix_sums(adata2) if adata2 is not None else {}

    # If only one dataset, use sums1
    label_to_use = ""
    if adata2 is None or not sums2:
        sums2 = {}
        label_to_use = "STARsolo"
    if adata1 is None or not sums1:
        sums1 = {}
        label_to_use = "kallisto"
    # Determine which suffixes to include (union)
    all_suffixes = sorted(set(sums1) | set(sums2),
                          key=lambda k: max(sums1.get(k,0), sums2.get(k,0)),
                          reverse=True)

    # Only show if 2–100 suffixes
    if not (2 <= len(all_suffixes) <= 100):
        return ""

    # Build traces
    if sums1:
        trace1 = go.Bar(
            x=[sums1.get(s,0) for s in all_suffixes],
            y=all_suffixes,
            orientation="h",
            name=label1,
            marker=dict(color="teal"),
            visible=True
        )
        data_traces = [trace1]

    updatemenus = []
    if sums2:
        trace2 = go.Bar(
            x=[sums2.get(s,0) for s in all_suffixes],
            y=all_suffixes,
            orientation="h",
            name=label2,
            marker=dict(color="teal"),
            visible=False if sums1 else True
        )
        if sums1:
            data_traces.append(trace2)
        else:
            data_traces = [trace2]

        if sums1:
        # Toggle buttons
            buttons = [
                {
                    "method": "update",
                    "label": label1,
                    "args": [{"visible": [True, False]}, {}]
                },
                {
                    "method": "update",
                    "label": label2,
                    "args": [{"visible": [False, True]}, {}]
                }
            ]
            updatemenus = [{
                "buttons": buttons,
                "direction": "right",
                "x": 0, "y": 1.015,
                "xanchor": "left", "yanchor": "top"
            }]

    # Layout
    layout = go.Layout(
        barmode="overlay",
        width=800, height=len(all_suffixes)*25 + 100,
        updatemenus=updatemenus,
        xaxis=dict(title=f"UMI count {label_to_use}"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)',
        margin=dict(t=20, l=50, r=20, b=50)
    )
    fig = go.Figure(data=data_traces, layout=layout)
    div = plot_div_once(fig)

    return div

