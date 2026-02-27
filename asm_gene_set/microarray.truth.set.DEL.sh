#!/bin/bash

#SBATCH -J truthsetDEL                 # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 1:00:00                      # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -n 8                            # Number of cores

# CB_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/truth_set_analysis/CB4856"
# CB_microarray="$CB_dir/CB4856_microarray_deletion.noHDRs.tsv"

# CB_files=(
#     $CB_dir/CB4856.paftools.DEL.tsv
#     $CB_dir/CB4856.pav.DEL.correctedEND.tsv
#     $CB_dir/CB4856.sniffles.DEL.tsv
# )

# echo -e 'CHROM\tPOS\tpaftools\tPAV\tsniffles' > $CB_dir/CB_header.tsv

# for SV_anno_file in ${CB_files[@]}; do
#     basename=$(basename $SV_anno_file)   
#     toolname=${basename#CB4856.}              
#     toolname=${toolname%.DEL*}   

#     output=$CB_dir/${toolname}.DEL.corroboration.tsv
    
#     awk -F'\t' '
#     NR==FNR { sv_start[$1, $2] = $2; sv_end[$1, $2] = $3; next } 
#     {
#         found = 0;
#         for (key in sv_start) {
#             split(key, coords, SUBSEP);
#             if ($1 == coords[1] && $2 >= sv_start[key] && $2 <= sv_end[key]) {
#                 found = 1;
#                 break;
#             }
#         }
#         print found
#     }' $SV_anno_file $CB_microarray > $output
# done

# paste $CB_microarray $CB_dir/paftools.DEL.corroboration.tsv \
#     $CB_dir/pav.DEL.corroboration.tsv \
#     $CB_dir/sniffles.DEL.corroboration.tsv > $CB_dir/CB_merged_tools.tsv

# cat $CB_dir/CB_header.tsv $CB_dir/CB_merged_tools.tsv > $CB_dir/CB_final_corroboration.tsv




JU_dir="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/truth_set_analysis/JU258"
JU_microarray="$JU_dir/JU258_microarray_deletion.noHDRs.tsv"

JU_files=(
    $JU_dir/JU258.paftools.DEL.tsv
    $JU_dir/JU258.pav.DEL.correctedEND.tsv
    $JU_dir/JU258.sniffles.DEL.tsv
)


echo -e 'CHROM\tPOS\tpaftools\tPAV\tsniffles' > $JU_dir/JU_header.tsv

for SV_anno_file in ${JU_files[@]}; do
    basename=$(basename $SV_anno_file)   
    toolname=${basename#JU258.}              
    toolname=${toolname%.DEL*}   

    output=$JU_dir/${toolname}.DEL.corroboration.tsv
    
    awk -F'\t' '
    NR==FNR { sv_start[$1, $2] = $2; sv_end[$1, $2] = $3; next } 
    {
        found = 0;
        for (key in sv_start) {
            split(key, coords, SUBSEP);
            if ($1 == coords[1] && $2 >= sv_start[key] && $2 <= sv_end[key]) {
                found = 1;
                break;
            }
        }
        print found
    }' $SV_anno_file $JU_microarray > $output
done

paste $JU_microarray $JU_dir/paftools.DEL.corroboration.tsv \
    $JU_dir/pav.DEL.corroboration.tsv \
    $JU_dir/sniffles.DEL.corroboration.tsv > $JU_dir/JU_merged_tools.tsv

cat $JU_dir/JU_header.tsv $JU_dir/JU_merged_tools.tsv > $JU_dir/JU_final_corroboration.tsv



## Removing deletion coordinates that fall within HDRs 
# zgrep "JU258" /vast/eande106/data/c_elegans/WI/divergent_regions/20231213/20231213_c_elegans_divergent_regions_strain.bed.gz > /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/truth_set_analysis/JU258/JU_HDRs.bed

# awk -F'\t' '
# NR==FNR { 
#     # Store start and end positions for each chromosome from the BED file
#     start[$1, NR] = $2; 
#     end[$1, NR] = $3; 
#     next 
# } 
# {
#     # Check if the position in CB_microarray falls within any BED region
#     exclude = 0;
#     for (i in start) {
#         split(i, key, SUBSEP);
#         if ($1 == key[1] && $2 >= start[i] && $2 <= end[i]) {
#             exclude = 1;
#             break;
#         }
#     }
#     if (!exclude) print $0;
# }' JU_HDRs.bed JU258_microarray_deletion.tsv > JU258_microarray_deletion.noHDRs.tsv

