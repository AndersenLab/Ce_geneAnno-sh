#!/bin/bash

#SBATCH -J helxier
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/helixerCB4856.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/helixerCB4856.rr 

module load singularity 

### TEST RUN ### 

# Download lineage files and example genome:
# singularity run $container_image fetch_helixer_models.py
    # retrieved list of available models from https://raw.githubusercontent.com/weberlab-hhu/Helixer/main/resources/model_list.csv
    # saved model land_plant_v0.3_a_0080.h5 to /home/loconn13/.local/share/Helixer/models
    # saved model vertebrate_v0.3_m_0080.h5 to /home/loconn13/.local/share/Helixer/models
    # saved model fungi_v0.3_a_0100.h5 to /home/loconn13/.local/share/Helixer/models
    # saved model invertebrate_v0.3_m_0100.h5 to /home/loconn13/.local/share/Helixer/models

# cd /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies/helixer_test
    # curl -L ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_lyrata/dna/Arabidopsis_lyrata.v.1.0.dna.chromosome.8.fa.gz --output Arabidopsis_lyrata.v.1.0.dna.chromosome.8.fa.gz

# output_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/gffs/helixer_test" 
# container_image="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images/loconn13999-helixer_default_v0.3.4.2025_01_30_updatedUser.sif"
# intput_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies/helixer_test"

# singularity run --bind $output_dir:/output_dir \
#     $container_image \
#     Helixer.py --fasta-path /output_dir/Arabidopsis_lyrata.v.1.0.dna.chromosome.8.fa.gz \
#     --lineage land_plant \
#     --species Arabidopsis_lyrata \
#     --gff-output-path /output_dir/Arabidopsis_lyrata_chromosome8_helixer.gff3

# ### C elegans run ###
output_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/gffs"
fasta_path="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies"
container_image="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/container_images/loconn13999-helixer_default_v0.3.4.2025_01_30_updatedUser.sif"

### GFFs to compare to
# /vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/c_elegans.PRJNA13758.WS283.annotations.gff3
# /vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.gff3
# /vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3
# /vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/cleanup_N2_WS283.sh - USE TO GENERATE FILTERED GFFS

singularity run --bind $output_dir:/output_dir --bind $fasta_path:/fasta_path \
    $container_image Helixer.py \
    --fasta-path /fasta_path/CB4856_hifi.inbred.asm.fa \
    --lineage invertebrate \
    --species Caenorhabditis_elegans \
    --gff-output-path /output_dir/CB4856.hifi.helixer.gff3