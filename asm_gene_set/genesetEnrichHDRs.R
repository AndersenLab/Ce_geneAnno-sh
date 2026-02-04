library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)
library(cowplot)

# ======================================================================================================================================================================================== #
# Loading in orthogroups
# ======================================================================================================================================================================================== #
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))

strainCol <- colnames(ortho_genes_dd)
ugh <- gsub(".20251012.inbred.blobFiltered.softMasked.braker.longestIso.protein","", strainCol)
ugh2 <- gsub(".20251014.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh)
ugh3 <- gsub(".20251124.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh2)
ugh4 <- gsub(".20251012.inbred.onlyONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh3)
ugh5 <- gsub(".Nov2025.softMasked.braker.longest.protein","", ugh4)
ugh6 <- gsub(".20251012.inbred.withONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh5)
strainCol_c2 <- gsub("c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein","N2", ugh6)
colnames(ortho_genes_dd) <- strainCol_c2

ortho_count <- ortho_genes_dd

strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations_pre <- ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))


private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups_UnassignedGenes.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))

colnames(private_OGs) <- strainCol_c2

private_cols <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup")]

private_ortho_count <- private_OGs
for (i in 1:length(private_cols)) {
  print(paste0(i, " out of ", length(private_cols)))
  temp_colname <- paste0(private_cols[i], "_count")
  
  private_ortho_count <- private_ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := ifelse(is.na(!!sym(private_cols[i])), NA, 1))
}

all_relations_private <- private_ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))

all_relations <- all_relations_pre %>%
  dplyr::bind_rows(all_relations_private)

private_freq = (1/(length(strainCol_c2_u)))

class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(Orthogroup,class)

all_class <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::left_join(class, by = "Orthogroup") 

long_class <- all_class %>%
  tidyr::pivot_longer(
    cols = -c(Orthogroup, class),
    names_to = "strain",
    values_to = "gene",
    values_drop_na = TRUE) %>%
  tidyr::separate_rows(gene, sep = ",\\s*") %>%
  dplyr::select(strain, gene, class) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene))

# test <- all_class %>% dplyr::group_by(class) %>% dplyr::mutate(count = n()) %>% dplyr::distinct(class,count)

# Loading all genes in pangenome
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/142_140WSs_andCGC1_longestIsoGenes.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
all_genes_strain <- rbind(genes_strain,N2_gff) %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::rename(gene = attributes)

genes_class <- all_genes_strain %>%
  dplyr::left_join(long_class, by = c("gene","strain"))


# HDRs
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
WSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/140_correct.tsv", col_names = "strain") %>% dplyr::pull()

# COLLAPSING HDRS AMONG 140 WSs
all_regions <- hdrs %>%
  dplyr::filter(strain %in% WSs) %>%
  dplyr::arrange(chrom,start) %>%
  data.table::as.data.table()


# ======================================================================================================================================================================================== #
# Extracting the longest WS contig alignment for every HDR #
# ======================================================================================================================================================================================== #
nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>% 
  dplyr::filter(strain != "ECA396") %>%
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) # add resolution for inverted alignments - need to pull genes differently

nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, WSS, WSE, contig, L2, LENQ, inv, strain) %>%
  data.table::as.data.table()

# Finding overlap of WS alignments with HDRs - for an INDIVIDUAL strain:
### Once finalalized, create a function with the analysis workflow and iterate over wild strains: for (i in WSs)

# Initialize a data frame and then add data to data frame for each strain
# columns: Strain, num_N2_hdrs, num_WS_aln_HDRs, mean_N2_hdr_size, mean_WS_hdr_size, max_WS_hdr_size, min_WS_hdr_size, num_WS_hdr_liftover
results_df = as.data.frame(matrix(ncol = 18, nrow = 140))
names(results_df) = c("Strain", "num_N2_hdrs", "num_WS_aln_HDRs", "mean_N2_hdr_size", "median_N2_hdr_size", 
                      "mean_WS_hdr_size", "median_WS_hdr_size", 'max_WS_hdr_size', "min_WS_hdr_size", "num_WS_hdr_liftover",
                      "num_WS_core_genes", "num_WS_acc_genes", "num_WS_priv_genes", "num_WS_core_inHDR", "num_WS_acc_inHDR", "num_WS_priv_inHDR", 
                      "total_WS_genes", "total_WS_genes_inHDRs")
for (i in 1:length(WSs)) {
  # print(i)
  SOI = WSs[i]
  # print(SOI)s
  results_df[i,1] = c(SOI)

  # Perform foverlaps
  
  
  
  
  results_df[i,2] = c(num_hdrs) # number of HDRs for strain i 
  results_df[i,3] = c(num_hdrs) # number HDRs that have WS alignment overlap
  
  # Wild strain N2 HDR stats
  results_df[i,4] = c() # mean N2 HDR size for strain i
  results_df[i,5] = c() # median N2 HDR size for strain i
  
  # Perform HDR liftover
  results_df[i,6] = c() # mean size of WS HDR liftover
  results_df[i,7] = c() # median size of WS HDR liftover
  results_df[i,8] = c() # maximum size of WS HDR liftover
  results_df[i,9] = c() # minumim size of WS HDR lifover
  results_df[i,10] = c() # number of WS HDR liftovers
  
  
  # Pull WS genes in lifted over HDRs
  
}















SOI <- "ECA3088"

nucmer_ranges <- nucmer_ranges %>% dplyr::filter(strain == SOI)
strain_hdr <- all_regions %>% dplyr::filter(strain == SOI) %>% 
  dplyr::rename(og_hdr_start = start, og_hdr_end = end) %>% 
  dplyr::mutate(start = ifelse(og_hdr_start >= 5000, og_hdr_start - 5000, og_hdr_start), end = og_hdr_end + 5000) ##################  EXTENDING HDR BOUNDARIES BY 5 KB ON BOTH SIDES TO TRY AND PULL EXTREMELY DIVERGENT REGIONS
  
strain_gene_class <- genes_class %>% dplyr::filter(strain == SOI) %>%
  dplyr::group_by(strain, class) %>%
  dplyr::mutate(strain_class_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(seqid,start,end,gene,strain,class,strain_class_count)

data.table::setkey(nucmer_ranges, chrom, start, end)
data.table::setkey(strain_hdr, chrom, start, end)

hdr_aln <- data.table::foverlaps(
  x = strain_hdr,
  y = nucmer_ranges,
  type = "any" # if the start/end of any alignment is within an HDR
  ) %>%
  dplyr::rename(hdr_start_extended = i.start, hdr_end_extended = i.end, N2S = start, N2E = end) %>%
  dplyr::select(-i.strain) # 10,333 genes!

num_hdrs <- nrow(strain_hdr)
num_hdr_aln <- nrow(hdr_aln %>% dplyr::distinct(hdr_start_extended))
strain_id <- hdr_aln %>% dplyr::distinct(strain) %>% dplyr::pull()
print(paste0(num_hdr_aln, " / ", num_hdrs, " HDRs have overlapping ", strain_id, " alignments"))

# Test plot
ggplot(data = hdr_aln) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

# ggplot(data = hdr_aln %>% dplyr::filter(chrom == "IV" & og_hdr_start > 12000000 & og_hdr_end < 14000000)) +
#   geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = contig), size = 1) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


nucmer_longest <- hdr_aln %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(og_hdr_start) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(og_hdr_start, contig) %>%
  dplyr::mutate(ntig = n()) %>%
  dplyr::mutate(tigsize = sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
  dplyr::ungroup() %>% 
  dplyr::group_by(og_hdr_start) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::group_by(og_hdr_start) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()

ggplot(data = nucmer_longest) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


# Example of farrrrrrr muli-alignments of the same contig:
jump <- nucmer_longest %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000")

ggplot(data = jump) +
  geom_rect(data = jump %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

# I need to remove very distanct alignments from the same contig
nucmer_longest_jumpRemoved <- nucmer_longest %>%
  dplyr::mutate(St2=ifelse(inv==T,WSE,WSS),Et2=ifelse(inv==T, WSS, WSE)) %>%
  dplyr::arrange(St2) %>%
  dplyr::group_by(chrom, og_hdr_start) %>%
  dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 1.5E5, 1 ,0)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom,og_hdr_start,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom,og_hdr_start) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup()

jump_rm <- nucmer_longest_jumpRemoved %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000")

ggplot(data = jump_rm) +
  geom_rect(data = jump_rm %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))



# Calculating the WS coordinates that directly overlap with the N2 HDR boundaries
nucmer_slope <- nucmer_longest_jumpRemoved %>%
  # dplyr::group_by(og_hdr_start) %>%
  # dplyr::mutate(longest_L2 = max(L2)) %>% # filtering for the longest alignment for each HDR to calculate the slope
  # dplyr::ungroup() %>%
  dplyr::mutate(slope = ((WSE - WSS) / (N2E - N2S))) %>%
  dplyr::mutate(intercept = WSS - (slope * N2S)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_hdr_start = ((slope * og_hdr_start) + intercept)) %>%
  dplyr::mutate(WS_hdr_end = ((slope * og_hdr_end) + intercept)) %>%
  dplyr::mutate(spans_hdr = (N2S <= og_hdr_start & N2E >= og_hdr_start) | (N2S <= og_hdr_end & N2E >= og_hdr_end)) %>%
  dplyr::group_by(chrom, og_hdr_start, spans_hdr) %>%
  # Add resolution on if there is only an alignment for one edge of the HDR
  dplyr::mutate(
    WS_hdr_start_min = ifelse(WSS == min(WSS) & inv == F & spans_hdr == T, WS_hdr_start, 
                              ifelse(WSE == min(WSE) & inv == T & spans_hdr == T, WS_hdr_end, NA))) %>% # conditional based on if alignment is INV and min WSS! - what if there is only an alignment for one HDR boundary???  # (N2S < og_hdr_start & N2E > og_hdr_start | N2S < og_hdr_end & N2E > og_hdr_end)
  dplyr::mutate(
    WS_hdr_end_max = ifelse(WSE == max(WSE) & inv == F & spans_hdr == T, WS_hdr_end, 
                            ifelse(WSS == max(WSS) & inv == T & spans_hdr == T, WS_hdr_start, NA))) %>% # conditional based on if alignment is INV and max WSE! - what if there is only an alignment for one HDR boundary???
  dplyr::ungroup() %>%
  dplyr::group_by(chrom, og_hdr_start) %>%
  tidyr::fill(WS_hdr_start_min, .direction = "updown") %>%
  tidyr::fill(WS_hdr_end_max, .direction = "updown") %>%
  dplyr::ungroup() 


large <- nucmer_slope %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000")

ggplot(data = large) +
  geom_rect(data = large %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = large %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  # coord_cartesian(ylim = c(0,5)) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

# After adding the condition that if an HDR has only one boundary with an alignment, and the other has an alignment withing X distance of the edge of the other alignment, then use the other alignment for calculating the slope
nucmer_slope_oneside <- nucmer_slope %>%
  dplyr::mutate(spans_hdr_start = (N2S <= og_hdr_start & N2E >= og_hdr_start),
                spans_hdr_end = (N2S <= og_hdr_end & N2E >= og_hdr_end)) %>%
  dplyr::mutate(outside_proximity_left = ifelse(spans_hdr_start == FALSE & N2E < og_hdr_start, og_hdr_start - N2E, NA)) %>%
  dplyr::mutate(outside_proximity_right = ifelse(spans_hdr_end == FALSE & N2S > og_hdr_end, N2S - og_hdr_end, NA)) %>%
  dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
  dplyr::mutate(WS_hdr_start_min_updated = ifelse(inv == T & WS_hdr_end != WS_hdr_start_min & !is.na(outside_proximity_right) & all(spans_hdr_end == F), WS_hdr_end, 
                                                  ifelse(inv == F & WS_hdr_start != WS_hdr_start_min & !is.na(outside_proximity_left) & all(spans_hdr_start == F), WS_hdr_start, NA))) %>%
  dplyr::mutate(WS_hdr_end_max_updated = ifelse(inv == T & WS_hdr_start != WS_hdr_end_max & !is.na(outside_proximity_left) & all(spans_hdr_start == F), WS_hdr_start, 
                                                  ifelse(inv == F & WS_hdr_end != WS_hdr_end_max & !is.na(outside_proximity_right) & all(spans_hdr_end == F), WS_hdr_end, NA))) %>%
  tidyr::fill(WS_hdr_start_min_updated, WS_hdr_end_max_updated, .direction = "downup") %>%
  dplyr::mutate(WS_hdr_start_min_updated = ifelse(is.na(WS_hdr_start_min_updated),WS_hdr_start_min,WS_hdr_start_min_updated),
                WS_hdr_end_max_updated = ifelse(is.na(WS_hdr_end_max_updated),WS_hdr_end_max,WS_hdr_end_max_updated)) %>%
  dplyr::ungroup()


# Account for situations where there are more than one alignment outside of a boundary for a single HDR!!
ggplot(data = nucmer_slope_oneside) + 
  geom_histogram(aes(x = outside_proximity_left)) +
  theme_classic()

ggplot(data = nucmer_slope_oneside) + 
  geom_histogram(aes(x = outside_proximity_right)) +
  theme_classic()



large_updated <- nucmer_slope_oneside %>% dplyr::filter(chrom == "IV" & og_hdr_start == "14170000") 

ggplot(data = large_updated) +
  geom_rect(data = large_updated %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = large_updated %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  # coord_cartesian(ylim = c(0,5)) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


# Now the other side!
test2 <- nucmer_slope_oneside %>% dplyr::filter(inv == F & !is.na(outside_proximity_right))

ugh <- nucmer_slope_oneside %>% dplyr::filter(chrom == "X" & og_hdr_start == "14355000")
ugh2 <- nucmer_slope %>% dplyr::filter(chrom == "X" & og_hdr_start == "14355000")

ggplot(data = ugh2) +
  geom_rect(data = ugh2 %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = ugh2 %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  # coord_cartesian(ylim = c(0,5)) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

ggplot(data = ugh) +
  geom_rect(data = ugh %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = ugh %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  # coord_cartesian(ylim = c(0,5)) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


# Plot HDRs that have no alignment spanning the edges!
no_aln_boundaries <- nucmer_slope_oneside %>% 
  dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
  dplyr::filter(all(spans_hdr == F))

ggplot(no_aln_boundaries %>% dplyr::filter(og_hdr_start == "19480000")) + 
  geom_rect(data = no_aln_boundaries %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = no_aln_boundaries %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

no_aln_bound_inv <- nucmer_slope_oneside %>% 
  dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
  dplyr::filter(all(spans_hdr == F)) %>%
  dplyr::filter(chrom == "V" & og_hdr_start == "17045000")

ggplot(no_aln_bound_inv) + 
  geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = no_aln_bound_inv %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))



nucmer_slope_final <- nucmer_slope_oneside %>%
  dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
  dplyr::mutate(WS_hdr_start_min_updated = ifelse(is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == F & WSS == min(WSS) & !is.na(outside_proximity_left) |
                                                    is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == F & WSE == max(WSE) & !is.na(outside_proximity_right), WS_hdr_start, 
                                                         ifelse(is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == T & WSE == min(WSE) & !is.na(outside_proximity_right) | 
                                                                  is.na(WS_hdr_start_min_updated) & all(spans_hdr == F) & inv == T & WSS == max(WSS) & !is.na(outside_proximity_left), WS_hdr_end, WS_hdr_start_min_updated))) %>%
  tidyr::fill(WS_hdr_start_min_updated, .direction = "updown") %>%
  dplyr::mutate(WS_hdr_end_max_updated = ifelse(is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == F & WSS == min(WSS) & !is.na(outside_proximity_left) | 
                                                  is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == F & WSE == max(WSE) & !is.na(outside_proximity_right), WS_hdr_end, 
                                                       ifelse(is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == T & WSE == min(WSE) & !is.na(outside_proximity_right) |
                                                                is.na(WS_hdr_end_max_updated) & all(spans_hdr == F) & inv == T & WSS == max(WSS) & !is.na(outside_proximity_left), WS_hdr_start, WS_hdr_end_max_updated))) %>%
  tidyr::fill(WS_hdr_end_max_updated, .direction = "updown") %>%
  dplyr::ungroup()



ggplot(nucmer_slope_final %>% dplyr::filter(og_hdr_start == "19480000")) + 
  geom_rect(data = nucmer_slope_final %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = nucmer_slope_final %>% dplyr::filter(og_hdr_start == "19480000") %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


ggplot(nucmer_slope_final %>%   dplyr::filter(chrom == "V" & og_hdr_start == "17045000")) + 
  geom_rect(data = nucmer_slope_final %>%   dplyr::filter(chrom == "V" & og_hdr_start == "17045000") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = nucmer_slope_final %>%   dplyr::filter(chrom == "V" & og_hdr_start == "17045000") %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
  

test <- nucmer_slope_final %>% dplyr::filter(chrom == "IV" & og_hdr_start > 12600000 & og_hdr_end < 13000000)
ggplot(test) + 
  geom_rect(data = test %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = test %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


test_inv <- nucmer_slope_final %>% dplyr::filter(chrom == "V" & og_hdr_start > 6200000 & og_hdr_end < 6700000)
ggplot(test_inv) + 
  geom_rect(data = test_inv %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = test_inv %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


ggplot(nucmer_slope_final) + 
  geom_rect(data = nucmer_slope_final %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = nucmer_slope_final %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


ggplot(nucmer_slope_final %>% dplyr::filter(chrom == "III" & N2S > 11000000)) +
  geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "III" & N2S > 11000000) %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = nucmer_slope_final %>% dplyr::filter(chrom == "III" & N2S > 11000000) %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


nucmer_slope_spanning_HDR_both_sides <- nucmer_slope_final %>%
  dplyr::filter(spans_hdr == T) %>%
  # dplyr::mutate(spans_hdr_start = (N2S <= og_hdr_start & N2E >= og_hdr_start),
                # spans_hdr_end = (N2S <= og_hdr_end & N2E >= og_hdr_end)) %>%
  dplyr::group_by(chrom,og_hdr_start, og_hdr_end) %>%
  dplyr::summarise(
    any_start = any(spans_hdr_start, na.rm = TRUE),
    any_end   = any(spans_hdr_end,   na.rm = TRUE), 
    n_align   = n(), .groups = "drop") %>%
  dplyr::filter(xor(any_start, any_end))

num_one_side_span <- nrow(nucmer_slope_spanning_HDR_both_sides)
print(paste0("Number of HDRs that only have an alignment overlapping with either the start or end HDR boundaries, but not both: ", num_one_side_span))


# nucmer_slope_spanning_HDR_both_sides_plt <- nucmer_slope_spanning_HDR_both_sides %>%
#   dplyr::left_join(nucmer_slope, by = c("chrom","og_hdr_start","og_hdr_end")) %>%
#   dplyr::group_by(og_hdr_start) %>%
#   dplyr::mutate(hdr_size = og_hdr_end - og_hdr_start, ws_hdr_size = WS_hdr_end_max - WS_hdr_start_min) 
# 
# 
# ggplot(data = nucmer_slope_spanning_HDR_both_sides_plt) +
#   geom_rect(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# ggplot(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::filter(chrom == "V")) +
#   geom_rect(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::filter(chrom == "V") %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::filter(chrom == "V") %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))
# 
# 
# ggplot(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::filter(chrom == "V" & N2S > 18000000 & N2E < 19500000)) +
#   geom_rect(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::filter(chrom == "V" & N2S > 18000000 & N2E < 19500000) %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = nucmer_slope_spanning_HDR_both_sides_plt %>% dplyr::filter(chrom == "V" & N2S > 18000000 & N2E < 19500000) %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


# only_one_boundary <- nucmer_slope_spanning_HDR_both_sides_plt %>%
#   group_by(chrom,og_hdr_start, og_hdr_end) %>%
#   summarise(ws_hdr_size = max(ws_hdr_size), .groups = "drop") %>%
#   slice_max(ws_hdr_size, n = 1)
# 
# only_one_boundary_plt <- nucmer_slope_spanning_HDR_both_sides_plt %>%
#   semi_join(only_one_boundary, by = c("chrom","og_hdr_start", "og_hdr_end"))
# 
# 
# ggplot(data = only_one_boundary_plt) +
#   geom_rect(data = only_one_boundary_plt %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = only_one_boundary_plt %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


# Finalizing set of WS HDR regionss
ws_n2_hdr_region <- nucmer_slope_final %>%
  dplyr::filter(spans_hdr == T) %>%
  dplyr::group_by(og_hdr_start) %>%
  dplyr::mutate(hdr_size = og_hdr_end - og_hdr_start, ws_hdr_size = WS_hdr_end_max_updated - WS_hdr_start_min_updated) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(og_hdr_start, .keep_all = T)

# Shows that HDRs are very structurally different than N2 and not just high SNV density
ggplot(data = ws_n2_hdr_region) +
  geom_point(aes(x = hdr_size / 1e6, y = abs(ws_hdr_size / 1e6)), color = 'firebrick') +
  geom_line(data = data.frame(x = c(0, 1.2)), aes(x = x, y = x), linetype = "dashed") +
  theme_bw() +
  labs(x = "N2 HDR size (Mb)", y = paste0(strain_id, " HDR size (Mb)"))

print(paste0("The largest HDR lift-over in ",SOI," is ", max(ws_n2_hdr_region$ws_hdr_size)))

# Distinct wild strain HDR lift-overs!
ws_hdr_coords <- nucmer_slope_final %>% 
  dplyr::distinct(chrom,N2S,N2E,og_hdr_start,og_hdr_end,longest_contig,WS_hdr_start_min_updated,WS_hdr_end_max_updated) %>%
  dplyr::filter(!is.na(WS_hdr_start_min_updated) & !is.na(WS_hdr_end_max_updated))

number_hdrs <- ws_hdr_coords %>% 
  dplyr::distinct(chrom,og_hdr_start,og_hdr_end, WS_hdr_start_min_updated, WS_hdr_end_max_updated) 


two_starts <- nucmer_slope_final %>% dplyr::filter(chrom == "II" & og_hdr_start == "1040000") %>%
  dplyr::select(chrom, og_hdr_start, og_hdr_end, longest_contig, N2S, N2E, WSS, WSE, WS_hdr_start_min_updated, WS_hdr_end_max_updated, spans_hdr, L2)

ggplot(data = two_starts) +
  geom_rect(data = two_starts %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = two_starts %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  # coord_cartesian(ylim = c(0,5)) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


aln_lengths <- nucmer_slope_final %>% dplyr::select(chrom, og_hdr_start, og_hdr_end, WS_hdr_start, WS_hdr_end, inv, WS_hdr_start_min_updated, WS_hdr_end_max_updated, L2) %>%
  dplyr::filter(!is.na(WS_hdr_start_min_updated) & !is.na(WS_hdr_end_max_updated)) %>%
  dplyr::filter((inv == T & WS_hdr_end == WS_hdr_start_min_updated | inv == T & WS_hdr_start == WS_hdr_end_max_updated) | 
                  (inv == F & WS_hdr_start == WS_hdr_start_min_updated | inv == F & WS_hdr_end == WS_hdr_end_max_updated))

number_hdrs_correct <- number_hdrs %>%
  dplyr::left_join(aln_lengths, by = c("chrom","og_hdr_start","og_hdr_end","WS_hdr_start_min_updated","WS_hdr_end_max_updated")) %>%
  dplyr::group_by(chrom, og_hdr_start) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom, og_hdr_start, og_hdr_end) %>%
  dplyr::filter(L2 == max(L2))

num <- nrow(number_hdrs_correct)

print(paste0("Total number of HDRs lifted over to wild strains: ", num, " / ", num_hdrs))

two_starts_after <- number_hdrs_correct %>% dplyr::filter(chrom == "II" & og_hdr_start == "1040000") %>%
  dplyr::left_join(two_starts, by = c("chrom","og_hdr_start","og_hdr_end", "WS_hdr_start_min_updated", "WS_hdr_end_max_updated"))
  # dplyr::filter(!is.na(L2.y))

ggplot(data = two_starts_after) +
  geom_rect(data = two_starts_after %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = two_starts_after %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  # coord_cartesian(ylim = c(0,5)) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


FINAL_hdrs <- number_hdrs_correct %>%
  dplyr::left_join(nucmer_slope_final, by = c("chrom","og_hdr_start","og_hdr_end", "WS_hdr_start_min_updated", "WS_hdr_end_max_updated")) %>%
  dplyr::mutate(WS_hdr_start_min_updated = ifelse(WS_hdr_start_min_updated < 0, 0, WS_hdr_start_min_updated))

# final_count <- FINAL_hdrs %>% dplyr::distinct(chrom,og_hdr_start)

# ughh <- ws_hdrs %>% dplyr::filter(chrom == "V" & og_hdr_start == "520000")
# ggplot(data = ughh) +
#   geom_rect(data = ughh %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(data = ughh %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
#   geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min_updated / 1e6, yend = WS_hdr_start_min_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max_updated / 1e6, yend = WS_hdr_end_max_updated / 1e6), color = 'black') +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = spans_hdr), size = 1) +
#   scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
#   # coord_cartesian(ylim = c(0,5)) +
#   # facet_wrap(~chrom) +
#   theme(
#     panel.border = element_rect(color = 'black', fill = NA),
#     panel.background = element_blank()
#   ) +
#   labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

ggplot(FINAL_hdrs) + 
  geom_rect(data = FINAL_hdrs %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = FINAL_hdrs %>% dplyr::distinct(WS_hdr_start_min_updated, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min_updated / 1e6, ymax = WS_hdr_end_max_updated / 1e6), fill = 'blue', alpha = 0.3) +
  # geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  # geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  # geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = spans_hdr), size = 1) +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "blue")) +
  facet_wrap(~chrom, scales = "free") +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


################## Look at HDRs that have "overlapping alignments" but then do not have HDR liftover ################################
final_count <- FINAL_hdrs %>% dplyr::distinct(chrom,og_hdr_start) # 313
pre_count <- nucmer_slope %>% dplyr::distinct(chrom,og_hdr_start) # 315

no_hdr_liftover <- dplyr::anti_join(nucmer_slope, FINAL_hdrs, by = c("chrom","og_hdr_start","og_hdr_end"))

ggplot(data = no_hdr_liftover) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


ggplot(data = nucmer_longest %>% dplyr::filter(og_hdr_start == "1269000")) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()) +
  # coord_cartesian(xlim = c(19.4,19.8)) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


ggplot(data = nucmer_longest %>% dplyr::filter(og_hdr_start == "13517000" & chrom == "II")) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()) +
  # coord_cartesian(xlim = c(19.4,19.8)) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))






# ======================================================================================================================================================================================== #
# Adding WS genes to each alignment, and then collapsing by N2 gene to get a matrix of syntenic genes #
# ======================================================================================================================================================================================== #


wsg <- data.table::as.data.table(ws_genes)
ws_hdrs <- FINAL_hdrs %>% dplyr::select(longest_contig, WS_hdr_start_min_updated, WS_hdr_end_max_updated, spans_hdr, chrom, og_hdr_start, og_hdr_end, N2S, N2E, WSS, WSE) %>%
  data.table::as.data.table()

# print() # the max and min WS start and end coordinate among any contig



clean_tigTrim <- tigTrim %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::select(chrom, unchanged_start_aln, unchanged_end_aln, start_aln, end_aln, L1, start_gene, end_gene, n2_gene, WSS, WSE, inv, longest_contig, nalign, L2, strain) %>%
  dplyr::mutate(invStartTrimmed = ifelse(inv == T, WSE, WSS), invEndTrimmed = ifelse(inv == T, WSS, WSE))

filt_nucm_long <- data.table::as.data.table(clean_tigTrim)

data.table::setnames(filt_nucm_long, c("longest_contig", "invStartTrimmed", "invEndTrimmed"), c("contig", "start", "end")) # use St2 Et2 because the gff coordinates will be reported in this way

data.table::setkey(wsg, contig, strain, start, end)
data.table::setkey(filt_nucm_long, contig, strain, start, end)

joined <- foverlaps(
  x = wsg,
  y = filt_nucm_long,
  type = "any" # if the start/end of the gene is contained within the N2 gene space or touching the boundaries of the WS alignment coordinates
) 

syntelog_matrix <- joined %>%
  dplyr::select(n2_gene,strain,attributes) %>%
  dplyr::rename(WS_gene = attributes) %>%
  tidyr::pivot_wider(
    id_cols = n2_gene,
    names_from = c(strain),
    values_from = c(WS_gene),
    values_fn = \(x) paste(unique(x), collapse = ",")
  ) %>%
  dplyr::filter(!is.na(n2_gene)) # remove the row that contains all non-syntenic predicted WS genes






































# Testing out Nic's code!
