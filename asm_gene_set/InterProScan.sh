#!/bin/bash

#SBATCH -J IntProScan                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 48                            # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/slurm.oe
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/slurm.rr

output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/output"
mkdir -p $output
TMP=${SLURM_TMPDIR:-/scratch4/eande106/Lance/ipr.$SLURM_JOB_ID}
mkdir -p $TMP

#         --applications Pfam,SMART,TIGRFAM,SUPERFAMILY,CDD,Gene3D,FunFam,PANTHER,PIRSF,PIRSR,ProSiteProfiles,ProSitePatterns,SFLD,Hamap,Coils,AntiFam \

/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/database/interproscan-5.75-106.0/interproscan.sh \
	--formats TSV \
	--input /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/prep/QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.protein.fa \
	--goterms \
	--cpu 48 \
	--iprlookup \
	--disable-precalc \
	--output-file-base $output/QX_IPR_allApps \
	--tempdir $TMP

