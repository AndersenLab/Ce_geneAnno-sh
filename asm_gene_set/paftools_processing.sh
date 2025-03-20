#!/bin/bash

#SBATCH -J paftoolsVCF                  # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 8:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 8                            # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/paf_manip.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/paf_manip.rr 

outdir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/updated_genomes"

for file in $outdir/*ALIGN.vcf; do
    strain=$(basename $file .hifi_asm20_1kbCOV_5kbALIGN.vcf)  
    bcftools sort $file -o ${file%.vcf}.sorted.vcf  
    bcftools reheader -s <(echo "$strain") -o ${file%.vcf}.sorted.renamed.vcf ${file%.vcf}.sorted.vcf
done

for file in $outdir/*sorted.renamed.vcf; do 
	bgzip -c $file > $file.gz
	tabix -p vcf $file.gz
done

bcftools merge -o $outdir/elegans.merged.1kbCOV.5kbALIGN.vcf -O v $outdir/*.sorted.renamed.vcf.gz

bcftools +fill-tags $outdir/elegans.merged.1kbCOV.5kbALIGN.vcf -- -t TYPE | bcftools view -o $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.vcf -O v

grep "^#CHROM" $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.vcf > $outdir/header.vcf
awk 'BEGIN {OFS="\t"} {print $1, $2, "end", $4, $5, $8, substr($0, index($0, $10))}' $outdir/header.vcf > $outdir/header_modified.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE[\t%GT]\n' $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.vcf | \
awk 'BEGIN { OFS="\t" } {
    if ($5 == "INDEL") {
        if (length($3) > length($4)) 
            $5 = "DEL";
        else if (length($3) < length($4)) 
            $5 = "INS";
    }
    print $0;
}' > $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.final0.vcf

awk 'BEGIN { OFS="\t" } {
    if ($5 == "SNP") 
        end = $2 + 1;
    else if ($5 == "INS") 
        end = $2 + length($4) - 1;
    else if ($5 == "DEL") 
        end = $2 + length($3) - 1;
    
    print $1, $2, end, $3, $4, $5, substr($0, index($0,$6));
}' $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.final0.vcf > $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.final1.vcf

cat $outdir/header_modified.vcf $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.final1.vcf > $outdir/elegans.merged.1kbCOV.5kbALIGN.annotated.final.vcf