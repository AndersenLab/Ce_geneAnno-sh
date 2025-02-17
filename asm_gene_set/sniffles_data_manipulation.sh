#!/bin/bash

#SBATCH -J snifflesdatamanipulation
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -c 6

### EXTRACTING INFORMATION FOR COMPARISON OF INDIVIDUAL VERSUS COMBINED MODE OF SNIFFLES2 
# grep -v "^#" combined35strain.sniffles2.vcf | awk -F'\t' '{print $1}' >> chrom.tsv
# grep -v "^#" combined35strain.sniffles2.vcf | awk -F'\t' '{print $2}' >> pos.tsv

# grep -v "^#" combined35strain.sniffles2.vcf | awk -F'\t' '
#     {
#         pos = $2;
#         end = pos + 1;  # Default to POS + 1 for BND

#         # Extract actual END if present
#         if ($8 ~ /END=[0-9]+/) {
#             match($8, /END=([0-9]+)/, arr);
#             end = arr[1];
#         }

#         # Adjust END for insertions
#         if ($8 ~ /SVTYPE=INS/ && $8 ~ /SVLEN=[-0-9]+/) {
#             match($8, /SVLEN=([-0-9]+)/, svlen_arr);
#             svlen = svlen_arr[1];
#             end = pos + svlen;  # Adjust END for insertions
#         }

#         print end;
#     }' >> end.tsv

# grep -v "^#" combined35strain.sniffles2.vcf | awk -F'\t' '{print $3}' | sed 's/Sniffles2.//' | sed 's/\..*//' >> annotation.tsv

# bcftools query -f '[%SAMPLE=%GT\t]\n' combined35strain.sniffles2.vcf | awk -F'\t' '{
#     ALT_samples = "";  
#     for (i = 1; i <= NF; i++) {  
#         if ($i ~ /0\/1|1\/0|1\/1/) { 
#             sub(/=.*/, "", $i);  
#             ALT_samples = ALT_samples (ALT_samples ? "," : "") $i; 
#         }
#     }
#     print (ALT_samples == "" ? "N/A" : ALT_samples);
# }' >> samples.tsv

# paste chrom.tsv pos.tsv end.tsv annotation.tsv samples.tsv > master_calls.tsv
# awk -F'\t' '$5 != "" {print $0}' master_calls.tsv >> master.tsv


## Now do for all individual VCFs!
for vcf in /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/sniffles2/vcf/*.vcf; do
    sample=$(basename $vcf .hifi_reads.vcf)  

    grep -v "^#" $vcf | awk -F'\t' '{print $1}' > ${sample}.chrom.tsv
   
    grep -v "^#" $vcf | awk -F'\t' '{print $2}' > ${sample}.pos.tsv

    grep -v "^#" $vcf | awk -F'\t' '
    {
        pos = $2;
        end = pos + 1;  # Default to POS + 1 for BND

        # Extract actual END if present
        if ($8 ~ /END=[0-9]+/) {
            match($8, /END=([0-9]+)/, arr);
            end = arr[1];
        }

        # Adjust END for insertions
        if ($8 ~ /SVTYPE=INS/ && $8 ~ /SVLEN=[-0-9]+/) {
            match($8, /SVLEN=([-0-9]+)/, svlen_arr);
            svlen = svlen_arr[1];
            end = pos + svlen;  # Adjust END for insertions
        }

        print end;
    }' >> ${sample}.end.tsv

    grep -v "^#" $vcf | awk -F'\t' '{print $3}' | sed 's/Sniffles2.//' | sed 's/\..*//' > ${sample}.annotation.tsv

    bcftools query -f '[%SAMPLE=%GT\n]' "$vcf" | awk -v strain="$sample" '
        {
            split($0, arr, "=");  # Split "SAMPLE=0/1" into arr[1]="SAMPLE" arr[2]="0/1"
            if (arr[2] == "0/1" || arr[2] == "1/1" || arr[2] == "1/0")
                print strain;
            else
                print "N/A";
        }' > ${sample}.ALT_samples.tsv

    paste ${sample}.chrom.tsv ${sample}.pos.tsv ${sample}.end.tsv ${sample}.annotation.tsv ${sample}.ALT_samples.tsv > /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/sniffles2/vcf/${sample}.SV.calls.tsv

    awk -F'\t' '$5 != "N/A" {print $0}' /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/sniffles2/vcf/${sample}.SV.calls.tsv >> /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/sniffles2/vcf/${sample}.SV.calls.final.tsv

    rm ${sample}.chrom.tsv ${sample}.pos.tsv ${sample}.end.tsv ${sample}.annotation.tsv ${sample}.ALT_samples.tsv

done