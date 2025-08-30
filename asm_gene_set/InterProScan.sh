#!/bin/bash

#SBATCH -J IntProScan                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 48                            # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/slurm.oe
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/slurm.rr


source activate interpro

output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/output"
mkdir -p $output
TMP=${SLURM_TMPDIR:-/scratch4/eande106/Lance/ipr.$SLURM_JOB_ID}
mkdir -p $TMP

interproscan.sh \
	--formats TSV \
	--input /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/prep/QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.protein.fa \
	--goterms \
	--cpu 48 \
	--applications Pfam,SMART,TIGRFAM,SUPERFAMILY,CDD,Gene3D \
	--iprlookup \
	--disable-precalc \
	--output-file-base $output/QX1410_InterProScan \
	--tempdir $TMP

