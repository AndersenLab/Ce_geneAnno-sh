#!/bin/bash

#SBATCH -J paftoolsSV
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/julAug2025_paftoolsRun.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/julAug2025_paftoolsRun.rr 
 
source activate 2024LabRetreat

# ref_genome_elegans="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
# asm_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies"
# out_paf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/paf/updated_genomes"
# out_vcf="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/updated_genomes"

# for hifi_asm in $asm_dir/*.hifi.inbred.fa; do
#     paf="$out_paf/$(basename $hifi_asm .inbred.fa)_asm20.paf"
#     vcf="$out_vcf/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
#     if [[ ! -f $vcf ]]; then
#         minimap2 -cx asm20 --cs -t 24 $ref_genome $hifi_asm > $paf # c flag - output paf
        
#         sort -k6,6 -k8,8n $paf \
#         | paftools.js call -l 1000 -L 5000 -f $ref_genome - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
#     fi
# done



# For 20205_07/08 briggsae
ref_genome_briggsae="/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa"
out_pafjulAug2025_CB="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/briggsae/paf/julAug2025"
out_vcfjulAug2025_CB="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/briggsae/vcf/julAug2025/"

while read -r species CB path_CB; do 
    paf="$out_pafjulAug2025_CB/$CB.hifi_asm20.paf"
    vcf="$out_vcfjulAug2025_CB/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 $ref_genome_briggsae $path_CB > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_briggsae - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/julaugCB.tsv"



# For 20205_07/08 tropicalis
ref_genome_tropicalis="/vast/eande106/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa"
out_pafjulAug2025_CT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/tropicalis/paf/julAug2025"
out_vcfjulAug2025_CT="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/tropicalis/vcf/julAug2025/"

while read -r species CT path_CT; do 
    paf="$out_pafjulAug2025_CT/$CT.hifi_asm20.paf"
    vcf="$out_vcfjulAug2025_CT/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 $ref_genome_tropicalis $path_CT > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_tropicalis - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/julaugCT.tsv"


# For 20205_07/08 elegans
ref_genome_elegans="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/c_elegans.PRJNA13758.WS283.genome.fa"
out_pafjulAug2025_CE="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/paf/julAug2025"
out_vcfjulAug2025_CE="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/vcf/julAug2025/"

while read -r species CE path_CE; do 
    paf="$out_pafjulAug2025_CE/$CE.hifi_asm20.paf"
    vcf="$out_vcfjulAug2025_CE/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 $ref_genome_elegans $path_CE > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_elegans - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/julaugCE.tsv"


# For 20205_07/08 nigoni
ref_genome_nigoni="/vast/eande106/data/c_nigoni/genomes/PRJNA384657/WS276/c_nigoni.PRJNA384657.WS276.genomic.fa.gz"
out_pafjulAug2025_CN="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/nigoni/paf/julAug2025"
out_vcfjulAug2025_CN="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/nigoni/vcf/julAug2025/"

while read -r species CN path_CN; do 
    paf="$out_pafjulAug2025_CN/$CN.hifi_asm20.paf"
    vcf="$out_vcfjulAug2025_CN/$(basename $paf .paf)_1kbCOV_5kbALIGN.vcf"
    
    if [[ ! -f $vcf ]]; then
        minimap2 -cx asm20 --cs -t 24 <(zcat $ref_genome_nigoni) $path_CN > $paf # c flag - output paf
        
        sort -k6,6 -k8,8n $paf \
        | paftools.js call -l 1000 -L 5000 -f $ref_genome_nigoni - > $vcf #-l min length to compute covereage, -L min alignment to call variants (default is 10kb and 50kb)
    fi
done < "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/julaugCN.tsv"









