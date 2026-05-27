#!/bin/bash

#SBATCH -J geneAnno                # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 12                           # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/matrixConcat.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/matrixConcat.rr 

#=================================

# Calling core, soft-core, rare, and private genes in Caenorhabditis elegans

#=================================

processed_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/gene_status"

# Concatenate all together to create a genotype matrix
key="$processed_data/1xCoverage/keyGene.txt"
awk -F',' '{print $1}' $processed_data/AB1.1x.98.genes.csv > $key

# mkdir -p $processed_data/1xCoverage 
# new="$processed_data/1xCoverage/fileGene.txt"
# awk -F',' '{print $1}' $processed_data/AB1.1x.98.genes.csv > $processed_data/1xCoverage/1xGenotypeMatrix.98.csv
# for file in $processed_data/*.1x.98.genes.csv; do
#     awk -F',' '{print $1}' $file > $new

#     if cmp -s $key $new; then
#         paste -d, $processed_data/1xCoverage/1xGenotypeMatrix.98.csv <(awk -F',' '{print $2}' $file) > temp_matrix.csv 
#         mv temp_matrix.csv $processed_data/1xCoverage/1xGenotypeMatrix.98.csv
#         rm $new
#     else   
#         echo "there is a difference in genes"
#     fi
# done

mkdir -p $processed_data/2xCoverage 
new="$processed_data/2xCoverage/fileGene.txt"
awk -F',' '{print $1}' $processed_data/AB1.2x.98.genes.csv > $processed_data/2xCoverage/2xGenotypeMatrix.98.csv
for file in $processed_data/*.2x.98.genes.csv; do
    awk -F',' '{print $1}' $file > $new

    if cmp -s $key $new; then
        paste -d, $processed_data/2xCoverage/2xGenotypeMatrix.98.csv <(awk -F',' '{print $2}' $file) > temp_matrix.csv 
        mv temp_matrix.csv $processed_data/2xCoverage/2xGenotypeMatrix.98.csv
        rm $new
    else   
        echo "there is a difference in genes"
    fi
done

mkdir -p $processed_data/3xCoverage 
new="$processed_data/3xCoverage/fileGene.txt"
awk -F',' '{print $1}' $processed_data/AB1.3x.98.genes.csv > $processed_data/3xCoverage/3xGenotypeMatrix.98.csv
for file in $processed_data/*.3x.98.genes.csv; do
    awk -F',' '{print $1}' $file > $new

    if cmp -s $key $new; then
        paste -d, $processed_data/3xCoverage/3xGenotypeMatrix.98.csv <(awk -F',' '{print $2}' $file) > temp_matrix.csv 
        mv temp_matrix.csv $processed_data/3xCoverage/3xGenotypeMatrix.98.csv
        rm $new
    else   
        echo "there is a difference in genes"
    fi
done

mkdir -p $processed_data/4xCoverage 
new="$processed_data/4xCoverage/fileGene.txt"
awk -F',' '{print $1}' $processed_data/AB1.4x.98.genes.csv > $processed_data/4xCoverage/4xGenotypeMatrix.98.csv
for file in $processed_data/*.4x.98.genes.csv; do
    awk -F',' '{print $1}' $file > $new

    if cmp -s $key $new; then
        paste -d, $processed_data/4xCoverage/4xGenotypeMatrix.98.csv <(awk -F',' '{print $2}' $file) > temp_matrix.csv 
        mv temp_matrix.csv $processed_data/4xCoverage/4xGenotypeMatrix.98.csv
        rm $new
    else   
        echo "there is a difference in genes"
    fi
done

mkdir -p $processed_data/5xCoverage 
new="$processed_data/5xCoverage/fileGene.txt"
awk -F',' '{print $1}' $processed_data/AB1.5x.98.genes.csv > $processed_data/5xCoverage/5xGenotypeMatrix.98.csv
for file in $processed_data/*.5x.98.genes.csv; do
    awk -F',' '{print $1}' $file > $new

    if cmp -s $key $new; then
        paste -d, $processed_data/5xCoverage/5xGenotypeMatrix.98.csv <(awk -F',' '{print $2}' $file) > temp_matrix.csv 
        mv temp_matrix.csv $processed_data/5xCoverage/5xGenotypeMatrix.98.csv
        rm $new
    else   
        echo "there is a difference in genes"
    fi
done

rm $key


