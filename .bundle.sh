#!/bin/sh

rm -rf swiftseq_pipeline
mkdir swiftseq_pipeline
cp config.txt swiftseq_pipeline/
# cp config_v2.txt swiftseq_pipeline/
# cp config_umirtbc.txt swiftseq_pipeline/
cp config.yaml swiftseq_pipeline/
cp cluster.yaml swiftseq_pipeline/
cp -R example_data swiftseq_pipeline/
cp example_samples.json swiftseq_pipeline/
cp -R scripts swiftseq_pipeline/
cp -R envs swiftseq_pipeline/
cp Snakefile swiftseq_pipeline/
cp run_pipeline.sh swiftseq_pipeline/
cp README.txt swiftseq_pipeline/
cp bundle.sh swiftseq_pipeline/.bundle.sh
cp sample_notebook.ipynb swiftseq_pipeline/
cp index_info.txt swiftseq_pipeline/

snakemake --forceall --dag | dot -Tpdf > dag.pdf
snakemake --forceall --rulegraph | dot -Tpdf > rulegraph.pdf

cp dag.pdf swiftseq_pipeline/
cp rulegraph.pdf swiftseq_pipeline/

date > swiftseq_pipeline/timestamp.txt

