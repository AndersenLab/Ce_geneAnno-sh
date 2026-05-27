#!/bin/bash

#SBATCH -J Ortho_Tran_gene                # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 24:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 8                           # Number of cores

sed -i 's|transcript_||g' /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_N2_braker_test.tsv

awk 'BEGIN { FS=OFS="\t" }
    NR==FNR { map[$1]=$2; next }  # Read key file: transcript -> gene
    {
        for (i=1; i<=NF; i++) {
            # Replace each comma-separated transcript with the gene(s)
	    n=split($i, a, /, */)
            for (j=1; j<=n; j++) {
                a[j] = (a[j] in map ? map[a[j]] : a[j])
            }
            $i = a[1]
            for (j=2; j<=n; j++) {
                $i = $i ", " a[j]
            }
        }
        print
    }' /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/transcripts_gene_115WI_N2.tsv /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_N2_braker_test.tsv > /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_N2_braker_genes.tsv

#awk -F'\t' 'BEGIN {OFS="\t"}
#NR == 1 {print; next}
#{
#  for (i = 4; i <= NF; i++) {
#    split($i, a, /, */)
#    delete seen
#    out = ""
#    for (j in a) {
#      if (a[j] != "" && !(a[j] in seen)) {
#        seen[a[j]] = 1
#        out = (out ? out ", " : "") a[j]
#      }
#    }
#    $i = out
#  }
#  print
#}' /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0425_genes.tsv > /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0425_genes_dedup.tsv
