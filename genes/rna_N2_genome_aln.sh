#!/bin/bash

#SBATCH -J RNA_aln
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=20G
#SBATCH --array=1-141%20

file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/rna_WS_genome_aln/temp_asm_list.tsv)

source activate star

base=$(basename $file)
strain=${base%%.*}

RNA_dir="/vast/eande106/data/c_elegans/WI/fastq/rna/eQTL_2019/merge_fq_2026"

if [[ ! -f $RNA_dir/$strain/${strain}.merged_r1.fq.gz ]]; then
    echo "RNA files not found for $strain, skipping"
    exit 0
fi

echo "Aligning RNA to genome for strain $strain ..."

OUT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/rna_WS_genome_aln"
mkdir -p $OUT/STAR_genome_index/$strain
mkdir -p $OUT/STAR_output/$strain

# Genome indexing
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 600000000000 \
--genomeDir $OUT/STAR_genome_index/$strain \
--genomeFastaFiles $file \
--genomeSAindexNbases 12

# Alignment
STAR \
--runThreadN 24 \
--genomeDir $OUT/STAR_genome_index/$strain \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--readFilesCommand zcat \
--alignIntronMax 10000 \
--limitBAMsortRAM 5000000000 \
--outSAMstrandField intronMotif \
--outFileNamePrefix $OUT/STAR_output/$strain/${strain}_ \
--readFilesIn $RNA_dir/$strain/${strain}.merged_r1.fq.gz $RNA_dir/$strain/${strain}.merged_r2.fq.gz

echo "Done with strain $strain"
