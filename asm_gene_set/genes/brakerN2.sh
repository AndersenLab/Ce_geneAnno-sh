#!/bin/bash

#SBATCH -J braker3
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24

module load singularity 

input="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/N2_BRAKER_error"
output="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/N2_BRAKER_error/container_braker3/output"
#busco="/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/busco_downloads/lineages"
container="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/software/containers/loconn13999-braker3_20250724.sif"
augustus_path="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/N2_BRAKER_error/container_braker3/augustus_config"

singularity exec \
	--bind $input:/input \
	--bind $output:/output \
	--bind $augustus_path:/augustus_config \
	--env AUGUSTUS_CONFIG_PATH=/augustus_config \
	$container \
	braker.pl \
	--genome /input/c_elegans.PRJNA13758.WS283.genome.fa \
	--species Ce_N2 \
	--prot_seq /input/N2_WS283.protein.fa \
	--threads 24 \
	--busco_lineage=nematoda_odb10 \
	--gff3 \
	--workingdir /output 
