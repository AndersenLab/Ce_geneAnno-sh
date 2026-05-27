#!/bin/bash

#SBATCH -J minimapLR
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/minimapLR.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/minimapLR.rr 
 
# for file in *.hifi_reads.bam; do sbatch --export=reads=$file /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/hifi_alignmentMinimap2.sh; done

source activate 2024LabRetreat

ref_genome="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
hifi_reads="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/lr_hifi"
out_aligned="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/lr_hifi/alignments"
out_fastq="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/lr_hifi/fastq"


strain=$(basename $reads .bam) 

samtools fastq -@ 24 $reads | gzip > $out_fastq/$strain.fq.gz

minimap2 -t 24 -ax map-hifi $ref_genome $out_fastq/$strain.fq.gz > $out_aligned/$strain.sam

samtools view -bS $out_aligned/$strain.sam | samtools sort -@ 24 -o $out_aligned/$strain.sorted.bam
samtools index $out_aligned/$strain.sorted.bam

rm $out_aligned/$strain.sam
