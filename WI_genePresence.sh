#!/bin/bash

#SBATCH -J geneAnno                # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 24                           # Number of cores

#=================================

# Calling core, soft-core, rare, and private genes in Caenorhabditis elegans

#=================================

# cd /vast/eande106/data/c_elegans/WI/alignments
# ls *.bam > /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/Ce_WI_allStrainList.txt

# while IFS= read -r strain; do sbatch --export=strain=${strain%.bam} getVariantCounts.sh; done < Ce_WI_allStrainList.txt 1 
# $1 is the coverage depth for a gene 

# mamba create -n gene_annotation mosdepth 
source activate gene_annotation

bams="/vast/eande106/data/c_elegans/WI/alignments"
raw_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data"
processed_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data"

ref_genome="c_elegans.PRJNA13758.WS283.genome.fa"
GFF="c_elegans.PRJNA13758.WS283.csq.gff3"

# strain_genome_cov_file="/data/eande106/eande106/analysis/alignment-nf/20231203-CE/20240213_c_elegans_gatk_ss.tsv" #1618 strains - 13 short.... there are 1631 total
# WIs_coverage=$processed_data/WIs_mean_genome_cov.csv

# Make 1kb windows of the reference genome
bedtools makewindows -g $raw_data/$ref_genome -w 1000 > $processed_data/N2.WS238.1kb_windows.bed

# Extract gene coordinates from the gff3
awk '$3 == "gene" {print $1, $4, $5, $9}' OFS="\t" $GFF | sed 's/ID=gene://;s/;.*//' > $processed_data/refGenesCoordinates.bed

# Intersect to get genome contents at each gene coordinate
bedtools intersect $processed_data/N2.WS238.1kb_windows.bed $processed_data/refGenesCoordinates.bed > $processed_data/N2.refGenes.bed

# Calculating coverage of each strain for ref genome
cd $processed_data/gene_coverage
mosdepth $strain $bams/$strain.bam --by $processed_data/N2.refGenes.bed --threads 24 --thresholds 1,2,3,4,5
# will output $strain.regions.bed.gz








cd $processed_data/gene_coverage
for WI_bam in $bams/*.bam; do
    mosdepth --by $processed_data/N2.refGenes.bed --fast-mode --threads 48 ${WI_bam%.bam} $WI_bam #output is ${WI_bam%bam}.regions.bed.gz. Also run without --fast-mode?? 
done

# CHANGE APPROACH TO BIN THE GENOME into 1kb bins and extract coverage > to a bed 
# then bedtools intersect of gene.ed with genome_coverage_1kbBIN.bed to get coverage at each gene
# then look at coverage of WI.bam for the bins and assess 1x, 2x, 3x, 4x, 5x for 80% or more of the bases in that gene and plot number of genes classified as present/absent and best coverage depth parameters is when the #s between two coverages are about the same - where the count stops drastically falling off
$1 is the coverage depth for a gene 

echo "STRAIN,AVE_COV" > $WIs_coverage
for coverage in $processed_data/gene_coverage/*.regions.bed.gz; do
    strain=${coverage%%.*} # Get strain name from filename
    strain_coverage=$(awk -v strain="$strain" '$1 == strain {print $4}' $strain_genome_cov_file) # Get average genome coverage for the strain

    echo "$strain,$strain_coverage" >> $WIs_coverage
    
    # Process gene coverage for this strain - ADD THE STRAIN NAME TO THE HEADER OF EACH FINAL FOR FINAL MERGING
    zcat $coverage | awk -v strain=$strain -v avg_cov=$strain_coverage '
        {
            # Calculate coverage ratio for each gene
            gene_cov_ratio = $4 / avg_cov

            if (gene_cov_ratio >= 0.8) {   # Criterion that the coverage for a gene must be at least 80% of the average genome-wide coverage for its strain
                print $1, $2, $3, $4, "1"  # Present 
            } else {
                print $1, $2, $3, $4, "0"  # Absent 
            }
        }' > $processed_data/gene_status/${strain}.gene_status.$1.tsv
done

final_matrix=$processed_data/final_matrix/gene_presence_matrix.$1.tsv 

# Aggregate results - each row is a gene, each column is a WI, and then each cell is a PRESENT (1) / ABSENT (0) for that gene in that strain
paste $processed_data/gene_status/*.gene_status.tsv | awk '{for (i=5; i<=NF; i+=5) printf "%s ", $i; print ""}' >> $final_matrix









### Write a python script to calculate frequency of a gene across all WIs
import pandas as pd
import matplotlib.pyplot as plt

# Load presence/absence matrix
df = pd.read_csv("gene_presence_matrix.tsv", sep="\t", header=None)

# Compute frequency of each gene across strains
freq = df.mean(axis=1) * 100  # Percentage of strains with the gene

# Classify genes
df['Category'] = pd.cut(freq, bins=[0, 5, 95, 100], labels=['rare', 'soft_core', 'core'])
df.loc[freq == 0, 'Category'] = 'Private'

df.to_csv("gene_classification.tsv", sep="\t", index=False)

plot the distrubition of frequency of genes appearing in WSs???? 
plot a pie chart representing the proportion of each type of gene???
