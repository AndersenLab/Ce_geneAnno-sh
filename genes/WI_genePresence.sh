#!/bin/bash

#SBATCH -J geneAnno                # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 12                           # Number of cores
#SBATCH --output=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/status.oe  
#SBATCH --error=/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/SLURM_output/status.rr 


#=================================

# Calling core, soft-core, rare, and private genes in Caenorhabditis elegans

#=================================

# cd /vast/eande106/data/c_elegans/WI/alignments
# ls *.bam > /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/Ce_WI_allStrainList.txt

# while IFS= read -r strain; do sbatch --export=strain=${strain%.bam} WI_genePresence.sh; done < ../../raw_data/Ce_WI_allStrainList.txt

# mamba create -n gene_annotation mosdepth 
source activate gene_annotation

bams="/vast/eande106/data/c_elegans/WI/alignments"
raw_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data"
processed_data="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data"

# cd $processed_data/gene_coverage
# mkdir -p $strain 
# cd $strain
# mosdepth $strain $bams/$strain.bam --by $processed_data/refGenesCoordinates.bed --threads 12 --thresholds 1,2,3,4,5

# Assessing presence/absence of genes
for coverage in 1 2 3 4 5; do

    echo "WBGeneID,$strain" > $processed_data/gene_status/$strain.${coverage}x.98.genes.csv

    zcat $processed_data/gene_coverage/$strain/$strain.thresholds.bed.gz | \
    awk -v coverage=$coverage -v threshold=0.98 '
    NR > 1 {
        gene = $4;                            
        total_bases = $3 - $2;                
        bases_covered = $(coverage + 4);     # Select the correct coverage column ($5 for 1x, $6 for 2x..... $9 for 5x)

        if ((bases_covered / total_bases) >= threshold) {
            print gene ",1";                 # Gene is present
        } else {
            print gene ",0";                 # Gene is absent
        }
    }' >> $processed_data/gene_status/$strain.${coverage}x.98.genes.csv
done


# # Concatenate all together to create a genotype matrix
# awk '{print $1}' 
# for file in *.1x.genes.csv; 
#     awk "$1 == $1 of XXX' then paste $2  
# done

# for file in *.2x.genes.csv; 
#     paste 
# done

# for file in *.3x.genes.csv; 
#     paste
# done 

# for file in *.4x.genes.csv; 
#     paste 
# done

# for file in *.5x.genes.csv; 
#     paste 
# done




# Calculating coverage of the binned reference genome
# cd $processed_data/binned_coverage
# mosdepth $strain $bams/$strain.bam --by $processed_data/N2.WS238.1kb_windows.bed --threads 12 --thresholds 1,2,3,4,5
# will output $strain.regions.bed.gz
# bedtools coverage -a $processed_data/refGenesCoordinates.bed -b $bams/$strain.bam -d > $processed_data/gene_coverage/$strain.refGenesCov.bedtoolsCoverage.bed
# bedtools coverage -a $processed_data/refGenesCoordinates.bed -b $bams/$strain.bam -hist > $processed_data/gene_coverage/$strain.refGenesCov.bedtoolsCoverageHist.bed





# # Intersect to get get the coverage at each gene locus
# zcat $processed_data/binned_coverage/$strain.regions.bed.gz | \
# bedtools intersect -a $processed_data/refGenesCoordinates.bed -b stdin -wa -wb > $processed_data/gene_coverage/$strain.refGenesCov.bed



# Reformatting and calculating fraction of gene that has coverage depth x
# awk 'BEGIN {
#     # Print headers
#     print "CHROM\tSTART\tEND\tGene\tTotal_Bases\tbasesCovered_1x\tFractionCovered_1x\tbasesCovered_2x\tFractionCovered_2x\tbasesCovered_3x\tFractionCovered_3x\tbasesCovered_4x\tFractionCovered_4x\tbasesCovered_5x\tFractionCovered_5x";
# } 
# {
#     gene = $4;
#     chrom[gene] = $1; # Store chromosome for the gene
#     start[gene] = $2; 
#     end[gene] = $3;
#     covered_1x[gene] += $9;
#     covered_2x[gene] += $10;
#     covered_3x[gene] += $11;
#     covered_4x[gene] += $12;
#     covered_5x[gene] += $13;
#     total_bases[gene] += ($3 - $2);
# } 
# END {
#     for (gene in total_bases) {
#         total = total_bases[gene];
#         frac_1x = (total > 0) ? covered_1x[gene] / total : 0;
#         frac_2x = (total > 0) ? covered_2x[gene] / total : 0;
#         frac_3x = (total > 0) ? covered_3x[gene] / total : 0;
#         frac_4x = (total > 0) ? covered_4x[gene] / total : 0;
#         frac_5x = (total > 0) ? covered_5x[gene] / total : 0;
#         print chrom[gene], start[gene], end[gene], gene, total, covered_1x[gene], frac_1x, covered_2x[gene], frac_2x, covered_3x[gene], frac_3x, covered_4x[gene], frac_4x, covered_5x[gene], frac_5x;
#     }
# }' OFS="\t" $processed_data/gene_coverage/$strain.refGenesCov.bedtoolsIntersect.bed > $processed_data/gene_coverage/$strain.CovSummary_refGenesBinnedGenome.tsv



# awk '
# BEGIN { OFS="\t"; print "CHROM", "START", "END", "WBGeneID", "total_bases", "1x_covered", "1x_fraction", "2x_covered", "2x_fraction", "3x_covered", "3x_fraction", "4x_covered", "4x_fraction", "5x_covered", "5x_fraction" }
# {
#     key = $1 FS $2 FS $3 FS $4
#     count[key]["total"]++
#     if ($6 >= 1) count[key]["1x"]++
#     if ($6 >= 2) count[key]["2x"]++
#     if ($6 >= 3) count[key]["3x"]++
#     if ($6 >= 4) count[key]["4x"]++
#     if ($6 >= 5) count[key]["5x"]++
# }
# END {
#     for (key in count) {
#         split(key, fields, FS)
#         total = count[key]["total"]
#         print fields[1], fields[2], fields[3], fields[4], total, \
#             count[key]["1x"], count[key]["1x"] / total, \
#             count[key]["2x"], count[key]["2x"] / total, \
#             count[key]["3x"], count[key]["3x"] / total, \
#             count[key]["4x"], count[key]["4x"] / total, \
#             count[key]["5x"], count[key]["5x"] / total
#     }
# }' gene_coverage.tsv > gene_thresholds_summary.tsv













