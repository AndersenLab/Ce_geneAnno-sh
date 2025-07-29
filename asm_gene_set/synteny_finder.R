library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(data.table)
library(cowplot)

# ======================================================================================================================================================================================== #
# Loading N0.tsv with transcripts converted to genes and plotting gene set classification #
# ======================================================================================================================================================================================== #
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/N0_0504_genes.tsv")

strainCol <- colnames(ortho_genes_dd)
strainCol_c1 <- gsub(".braker.longest.protein","",strainCol)
strainCol_c2 <- gsub(".longest.protein","",strainCol_c1)
colnames(ortho_genes_dd) <- strainCol_c2

ortho_count <- ortho_genes_dd

strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)

}
all_relations_pre <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))

private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/private_OG_genes_gene_assignments.tsv") %>%
  tidyr::pivot_wider(id_cols = Orthogroup, names_from = Source, values_from = Gene)

cols_yuh <- colnames(private_OGs)
cols_yuh1 <- gsub(".braker.longest.protein","", cols_yuh)
cols_yuh2 <- gsub(".longest.protein","", cols_yuh1)
colnames(private_OGs) <- cols_yuh2

private_cols <- cols_yuh2[!cols_yuh2 %in% c("Orthogroup")]

private_ortho_count <- private_OGs
for (i in 1:length(private_cols)) {
  print(paste0(i, " out of ", length(private_cols)))
  temp_colname <- paste0(private_cols[i], "_count")
  
  private_ortho_count <- private_ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := ifelse(is.na(!!sym(private_cols[i])), NA, 1))
}

all_relations_private <- private_ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count")) %>%
  dplyr::rename(HOG = Orthogroup)

all_relations <- all_relations_pre %>%
  dplyr::bind_rows(all_relations_private)

all_ortho_genes_dd <- ortho_genes_dd %>%
  dplyr::select(-OG, -'Gene Tree Parent Clade') %>%
  dplyr::bind_rows((private_OGs %>% dplyr::rename(HOG = Orthogroup)))

# Extracting HOG/OG identifiers for private OGs
private_freq = (1/(length(strainCol_c2_u)))

private_rowids <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>%
  dplyr::filter(class == "private") %>%
  dplyr::pull(HOG)


# ======================================================================================================================================================================================== #
# Splitting complex HOGs from 1-to-1s #
# ======================================================================================================================================================================================== #
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
all_genes_strain <- genes_strain %>%
  dplyr::filter(strain != "N2" | grepl("protein_coding", attributes)) %>% 
  dplyr::mutate(attributes = gsub("ID=gene:","", attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","", attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) 

ws_genes <- all_genes_strain %>% dplyr::filter(strain != "N2") %>% dplyr::select(-type)

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/nucmer_runs/115_WI_transformed_coords_FIXED.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY,-LENR)

all_relations_rowid <- all_relations %>% dplyr::rename(rowid = HOG)
ortho_genes_dd_rowid <- all_ortho_genes_dd %>% dplyr::rename(rowid = HOG)

# Extract rowids for simple and complex HOGs
# simple_rowids <- all_relations_rowid %>% 
#   dplyr::filter(if_all(2:ncol(.), ~ is.na(.) | . <= 1)) %>% 
#   dplyr::pull(rowid)

complex_all_rowids <- all_relations_rowid %>% 
  dplyr::filter(if_any(2:ncol(.), ~ . > 1)) %>% 
  dplyr::pull(rowid) 

complex_rowids <- setdiff(complex_all_rowids, private_rowids) #remove complex HOGs that are PRIVATE!

# simple_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% simple_rowids)
complex_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% complex_rowids) ### further decompress into only complex orthogroups that have ONE N2 gene

# noN2_complex_HOG <- complex_HOGS %>% dplyr::filter(is.na(N2))

# Creating a subset dataframe for BWC
# complex_rowids <- all_relations_rowid %>%
#   dplyr::mutate(row_sum = rowSums(dplyr::across(2:ncol(.)), na.rm = TRUE)) %>%
#   dplyr::filter(if_any(2:ncol(.), ~ . > 1)) %>%
#   dplyr::arrange(dplyr::desc(row_sum)) %>%
#   dplyr::slice_head(n = 10) %>%
#   dplyr::pull(rowid)
# 
# complex_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% complex_rowids) ### further decompress into only complex orthogroups that have ONE N2 gene
# 
# vis <- complex_HOGS %>%
#   dplyr::rename(HOG = rowid) %>%
#   dplyr::select(HOG,ECA1413,JU2617,XZ1515) %>%
#   dplyr::rename(Orthogroups = HOG)

n2_goi <- complex_HOGS%>%
  dplyr::select(N2) %>%
  tidyr::separate_rows(N2, sep = ",\\s*") %>%
  dplyr::filter(!is.na(N2)) 

n2_genes <- all_genes_strain %>% 
  dplyr::filter(strain == "N2") %>% 
  dplyr::rename(chrom = contig) %>%
  dplyr::filter(str_detect(attributes, str_c(n2_goi$N2, collapse = "|"))) %>%
  # dplyr::filter(reduce(n2_goi$N2, ~ .x | str_detect(attributes, fixed(.y)), .init = FALSE)) %>%
  data.table::as.data.table() # 10,333

# ======================================================================================================================================================================================== #
# Extracting the longest WS contig alignment for every N2 gene coordinate #
# ======================================================================================================================================================================================== #
nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, contig, WSS, WSE, L2, LENQ, strain) %>%
  data.table::as.data.table()

data.table::setkey(nucmer_ranges, chrom, start, end)
data.table::setkey(n2_genes, chrom, start, end)

n2_genes_align <- data.table::foverlaps(
  x = n2_genes,
  y = nucmer_ranges,
  type = "any" # if the start/end of any alignment is within an N2 gene
) %>%
  dplyr::rename(N2 = i.strain, n2_gene = attributes, start_aln = start, end_aln = end, start_gene = i.start, end_gene = i.end) %>%
  dplyr::select(-type, -strand) # 10,333 genes!


# Join where an N2 gene is fully contained within an alignment
# n2_genes_align <- fuzzyjoin::interval_inner_join(
#   x = n2_genes,
#   y = nucmer_ranges,
#   by = c("start", "end"),
#   type = "within"
# ) %>%
#   dplyr::filter(chrom.x == chrom.y) %>%
#   dplyr::rename(start_gene = start.x, end_gene = end.x, start_aln = start.y, end_aln = end.y, n2_gene = attributes, strain = strain.y, N2 = strain.x)
# # 10,284 N2 genes

nucmer_longest <- n2_genes_align %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, contig, strain) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
  dplyr::ungroup() %>% 
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()



# manyAln <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00008623")) + 
#   facet_wrap(~strain, scales = 'free', strip.position = "right") +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.05) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#         legend.position = 'none',
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_rect(fill = NA),
#         plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00008623 : IV")
# manyAln

# #### Testing how selection of "best" contig works with dot plots ####
# all_ctg_align <- n2_genes_align %>%
#   dplyr::group_by(strain, n2_gene) %>%
#   dplyr::mutate(nalign = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(n2_gene,contig, strain) %>%
#   dplyr::mutate(ntig= n()) %>%
#   dplyr::mutate(tigsize=sum(L2)) %>%
#   dplyr::ungroup()
# 
# all_test <- ggplot(nucmer_longest %>% dplyr::filter(strain == "ECA36" & n2_gene == "WBGene00010415")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.05) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA36 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00010415 : CHOROM")
# all_test
                     

# all_ctg_gene <- ggplot(all_ctg_align %>% dplyr::filter(n2_gene == "WBGene00013172" & strain == "ECA923")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.05) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA923 coord", x = "N2 coord", title = "WBGene00013172 : I")
# all_ctg_gene
# 
# lg_ctg_gene <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00013172" & strain == "ECA923")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.05) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA923 coord", x = "N2 coord", title = "WBGene00013172 : I")
# lg_ctg_gene
# 
# 
# all_ct_gene_test2 <- ggplot(all_ctg_align %>% dplyr::filter(n2_gene == "WBGene00016483" & strain == "QX1791")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "QX1791 coord", x = "N2 coord", title = "WBGene00016483 : V")
# all_ct_gene_test2
# 
# lg_ct_gene_test2 <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00016483" & strain == "QX1791")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "QX1791 coord", x = "N2 coord", title = "WBGene00016483 : V")
# lg_ct_gene_test2
# 
# 
# all_ct_gene_test3 <- ggplot(all_ctg_align %>% dplyr::filter(n2_gene == "WBGene00022291" & strain == "JU1581")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSE) / 1e6, ymax = max(WSS) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "JU1581 coord", x = "N2 coord", title = "WBGene00022291 : X")
# all_ct_gene_test3
# 
# lg_ct_gene_test3 <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00022291" & strain == "JU1581")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSE) / 1e6, ymax = max(WSS) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "JU1581 coord", x = "N2 coord", title = "WBGene00022291 : X")
# lg_ct_gene_test3
# 
# 
# all_ct_gene_test4 <- ggplot(all_ctg_align %>% dplyr::filter(n2_gene == "WBGene00271804" & strain == "JU311")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSE) / 1e6, ymax = max(WSS) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "JU1581 coord", x = "N2 coord", title = "WBGene00271804 : V")
# all_ct_gene_test4
# 
# lg_ct_gene_test4 <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00271804" & strain == "JU311")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSE) / 1e6, ymax = max(WSS) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "JU311 coord", x = "N2 coord", title = "WBGene00271804 : V")
# lg_ct_gene_test4
# 
# 
# all_ct_gene_test5 <- ggplot(all_ctg_align %>% dplyr::filter(n2_gene == "WBGene00269422" & strain == "ECA722")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSE) / 1e6, ymax = max(WSS) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA722 coord", x = "N2 coord", title = "WBGene00269422 : X")
# all_ct_gene_test5
# 
# lg_ct_gene_test5 <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00269422" & strain == "ECA722")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSE) / 1e6, ymax = max(WSS) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA722 coord", x = "N2 coord", title = "WBGene00269422 : X")
# lg_ct_gene_test5

###### Visualizing jumps for non single exon genes (>1kb) ###### 
# onekb <- nucmer_longest %>%
#   dplyr::filter((end_gene - start_gene) > 1000) # filtering for multi-exonic genes (>1kb)
# 
# g1plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00013300" & strain == "ECA369")) +
#     geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#     geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#     theme_bw() +
#     labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# g1plt
# 
# g2plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00018921" & strain == "MY10")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "MY10 coord", x = "N2 coord", title = "WBGene00018921 : I")
# g2plt
# 
# g3plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00018921" & strain == "ECA2676")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA2676 coord", x = "N2 coord", title = "WBGene00018921 : I")
# g3plt
# 
# g4plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00011254" & strain == "ECA1260")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA1260 coord", x = "N2 coord", title = "WBGene00011254 : V")
# g4plt
# 
# g5plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00004159" & strain == "NIC195")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "NIC195 coord", x = "N2 coord", title = "WBGene00004159 : IV")
# g5plt
#
# 
# g6plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00004159")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   facet_wrap(~strain, scales = "free", strip.position = "right") +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00004159")
# g6plt
# 
# 
# new <- nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00011254")
# g7plt <- ggplot(new) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.05) +
#   facet_wrap(~strain, scales = "free", strip.position = "right") +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   theme(
#         legend.position = 'none',
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_rect(fill = NA),
#         plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00011254 : V")
# g7plt

# 
# 
# g9plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00044801" & strain == "ECA1997")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA1997 coord", x = "N2 coord", title = "WBGene00044801 : IV")
# g9plt
# 
# 
# g10plt <- ggplot(nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00044801")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   facet_wrap(~strain, scales = "free", strip.position = "right") +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00044801 : III")
# g10plt



# ======================================================================================================================================================================================== #
# Removing contigs that are far away in WI coordinate system - large jumps in alignment #
# ======================================================================================================================================================================================== #
nucmer_longest_jump <- nucmer_longest %>%
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) %>% # addresses inverted alignments for tigTrim
  dplyr::arrange(St2) %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(leadDiff = lead(St2)-Et2) %>%
  dplyr::ungroup()
                
# ggplot(nucmer_longest_jump) +
#   geom_histogram(aes(x = leadDiff), bins = 300) +
#   geom_vline(xintercept = 4.5E5) +
#   geom_vline(xintercept = -4.5E5) +
#   # geom_vline(xintercept = -1.5E5) +
#   # coord_cartesian(x = c(-300000, 300000)) +
#   #scale_y_log10() +
#   theme_bw() +
#   xlim(-6E5,6E5)


# jmp <- ggplot(nucmer_longest_jumps %>% dplyr::filter(n2_gene == "WBGene00000838") %>% dplyr::filter(strain == "ECA2417")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   # geom_text(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00000838") %>% dplyr::filter(strain == "ECA2417") %>% dplyr::filter(leadDiff == max(leadDiff, na.rm = TRUE)),
#   # aes(x =  end_aln / 1e6 , y = Et2 / 1e6, label = paste0("Jump: ", leadDiff)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA2417 coord", x = "N2 coord", title = "WBGene00000838 : V")
# jmp


# jumpPlt <- ggplot(nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00020754") %>% dplyr::filter(strain == "ECA2529")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   geom_text(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00020754") %>% dplyr::filter(strain == "ECA2529") %>% dplyr::filter(leadDiff == max(leadDiff, na.rm = TRUE)),
#             aes(x =  end_aln / 1e6 , y = Et2 / 1e6, label = paste0("Jump: ", leadDiff)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#         legend.position = 'none',
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_rect(fill = NA),
#         plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA2529 coord", x = "N2 coord", title = "WBGene00020754 : V")
# jumpPlt
# 
# jumpPlt2 <- ggplot(nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00012264") %>% dplyr::filter(strain == "NIC199")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   geom_text(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00012264") %>% dplyr::filter(strain == "NIC199") %>% dplyr::filter(leadDiff == max(leadDiff, na.rm = TRUE)),
#             aes(x =  end_aln / 1e6 , y = Et2 / 1e6, label = paste0("Jump: ", leadDiff)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "NIC199 coord", x = "N2 coord", title = "WBGene00012264 : I")
# jumpPlt2
# 
# 
# jumpPlt_small <- ggplot(nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00010377") %>% dplyr::filter(strain == "XZ1513")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   geom_text(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00010377") %>% dplyr::filter(strain == "XZ1513") %>% dplyr::filter(leadDiff == max(leadDiff, na.rm = TRUE)),
#             aes(x =  end_aln / 1e6 , y = Et2 / 1e6, label = paste0("Jump: ", leadDiff)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "XZ1513 coord", x = "N2 coord", title = "WBGene00010377 : V")
# jumpPlt_small
# 
# jumpPlt_small2 <- ggplot(nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00004099") %>% dplyr::filter(strain == "JU775")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   geom_text(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00004099") %>% dplyr::filter(strain == "JU775") %>% dplyr::filter(leadDiff == max(leadDiff, na.rm = TRUE)),
#             aes(x =  end_aln / 1e6 , y = Et2 / 1e6, label = paste0("Jump: ", leadDiff)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "JU775 coord", x = "N2 coord", title = "WBGene00004099 : V")
# jumpPlt_small2

nucmer_longest_jumpRemoved <- nucmer_longest_jump %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(leadDiff = ifelse(is.na(leadDiff), 0, leadDiff)) %>%
  dplyr::mutate(jump = ifelse(abs(leadDiff) > 4.5E5, 1, 0)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, strain, run_id) %>%
  dplyr::mutate(gsize = n()) %>%
  dplyr::mutate(len = abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(sumlen == max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup()

#### Testing how my two-variable heuristics look when I do not remove the initial jump filter
nucmer_longest_jumpRemoved <- nucmer_longest_jump
#####

# check1 <- ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00013300" & strain == "ECA369")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = min(WSS) / 1e6, ymax = max(WSE) / 1e6), fill = "darkolivegreen4", alpha = 0.3) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# check1

# g10plt_jmRem <- ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00044801")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   facet_wrap(~strain, scales = "free", strip.position = "right") +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00044801")
# g10plt_jmRem
# 
# rmJumpPlt <- ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00020754") %>% dplyr::filter(strain == "ECA2529")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA2529 coord", x = "N2 coord", title = "WBGene00020754 : V")
# rmJumpPlt
# 
# rmJumpPlt2 <- ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00013300") %>% dplyr::filter(strain == "ECA369")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# rmJumpPlt2


### JUMPS AROUND 450kb ###
# fourfitty <- nucmer_longest_jumpRemoved %>%
#   dplyr::filter(jump == 0) %>%
#   dplyr::group_by(n2_gene, strain, L2) %>%
#   dplyr::filter(n() > 1) %>% 
#   dplyr::ungroup() %>%
#   dplyr::filter(leadDiff >= 0)
# 
# 
# 
# test_jmp_bound <- ggplot(fourfitty %>% dplyr::filter(n2_gene == "WBGene00001235") %>% dplyr::filter(strain == "ECA722")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   geom_text(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00001235") %>% dplyr::filter(strain == "ECA722") %>% dplyr::filter(leadDiff == max(leadDiff, na.rm = TRUE)),
#             aes(x =  end_aln / 1e6 , y = Et2 / 1e6, label = paste0("Jump: ", leadDiff)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA722 coord", x = "N2 coord", title = "WBGene00001235 : III")
# test_jmp_bound



## What is the distribution in contig lengths that are retained after filtering for long jumps? ##
# dist <- nucmer_longest_jumpRemoved %>%
#   dplyr::filter(n2_gene == "WBGene00044801")

# ctg_len_dist <- ggplot(dist) + 
#   # geom_bar(aes(x = L2), color = 'magenta4') +
#   geom_histogram(aes(x = L2), fill = 'magenta4', bins = 30) +
#   facet_wrap(~strain, scales = "free_x", strip.position = "top") +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA))
# ctg_len_dist
#   
# longest_ctg_dist <- dist %>%
#   dplyr::group_by(strain) %>%
#   dplyr::filter(L2 == max(L2)) %>%
#   dplyr::ungroup()
# 
# ctg_len_dist2 <- ggplot(lng_ctg %>% dplyr::filter(n2_gene == "WBGene00044801")) + 
#   # geom_bar(aes(x = L2), color = 'magenta4') +
#   geom_histogram(aes(x = L2), fill = 'magenta4', bins = 30) +
#   facet_wrap(~strain, scales = "free_x", strip.position = "top") +
#   theme(
#     legend.position = 'none',
#     axis.text.x = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA))
# ctg_len_dist2

# ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00011875")) + 
#   facet_wrap(~strain, scales = 'free', strip.position = "right") +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.05) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00011875 : III")
# 
# ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00011875" & strain == "ECA1208")) + 
#   # facet_wrap(~strain, Cscales = 'free', strip.position = "right") +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.05) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     # axis.text = element_blank(),
#     # axis.ticks = element_blank(),
#     # axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00011875 : III")


# ======================================================================================================================================================================================== #
# Selecting longest alignment for contigs that align twice after filtering for jumps 
# ======================================================================================================================================================================================== #
nucmer_longest_jumpRemoved <- nucmer_longest_jumpRemoved %>%
  dplyr::mutate(row_id = dplyr::row_number())

diff <- nucmer_longest_jumpRemoved %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(dplyr::n() == 2) %>%
  dplyr::arrange(L2, .by_group = TRUE) %>%
  dplyr::mutate(diff = abs(diff(L2)[1])) %>% 
  dplyr::mutate(diff_fraction = 1 - (min(L2) / max(L2))) %>%
  dplyr::arrange(St2, .by_group = TRUE) %>%
  dplyr::mutate(jumptwo = abs(lead(St2) - Et2)) %>%
  dplyr::ungroup() ### ~55,000 rows, so 27,500 different gene-strain entries, and then 237 genes per strain, so ~1 % of n2 genes for each strain

diff_dist_plot <- ggplot(data = diff) +
  geom_histogram(aes(x = diff / 1e3), bins = 500) + 
  theme_bw() +
  xlab("Pairwise difference in alignment (kb)")
diff_dist_plot       

diff_dist_plotjumps <- ggplot(data = diff) +
  geom_histogram(aes(x = jumptwo / 1e6), bins = 200, fill = 'black') +
  geom_vline(xintercept = max(diff$jumptwo, na.rm = TRUE) / 1e6, color = 'red', size = 2) +
  geom_vline(xintercept = 4.5, color = 'blue', size = 2) +
  # geom_vline(xintercept = 0.5, color = 'blue', size = 2) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  xlab("Difference in alignment jumps (absolute value) (Mb)")
diff_dist_plotjumps



# try_k <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00000864") %>% dplyr::filter(strain == "ECA768")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA768 coord", x = "N2 coord", title = "WBGene00000864 : II")
# try_k

# negJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00009423") %>% dplyr::filter(strain == "JU310")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "JU310 coord", x = "N2 coord", title = "WBGene00009423 : V")
# negJump
# 
# 
# zeroJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00000680") %>% dplyr::filter(strain == "PX179")) +
#     geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#     geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#     theme(
#       legend.position = 'none',
#       panel.background = element_blank(),
#       panel.grid = element_blank(),
#       panel.border = element_rect(fill = NA),
#       plot.title = element_text(size = 14, color = 'black')) +
#     labs(y = "PX179 coord", x = "N2 coord", title = "WBGene00000680 : I")
# zeroJump
# 
# lrgJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00000638") %>% dplyr::filter(strain == "XZ1513")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "XZ1513 coord", x = "N2 coord", title = "WBGene00000638 : I")
# lrgJump
# 
# # TWO DIFFERENT ALIGNMENTS FOR THE SAME CONTIG THAT ARE THE SAME SIZE
# largestJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00012066") %>% dplyr::filter(strain == "CB4856")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "CB4856 coord", x = "N2 coord", title = "WBGene00018439 : V")
# largestJump
#   
# contigDup <- ggplot() +
#   geom_segment(data = nucmer %>% dplyr::filter(contig == "ptg000004l") %>% dplyr::filter(strain == "CB4856"), aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   geom_rect(data = diff %>% dplyr::filter(strain == "CB4856" & n2_gene == "WBGene00012066"), aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "green") +
#   coord_cartesian(xlim = c(19.2,19.45)) + 
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "CB4856 coord", x = "N2 coord", title = "CB4856 ptg000004l")
# contigDup

# Removing alignments that are VERY small and distant from "main" alignment
# Need to find the WS coordinates that overlap with the N2 gene for each alignment, then calculate the mean (middle) WS coordinate, 
# and then if the difference is >trim_spacer (need to account for scale_distortion between WS and N2) and the smaller alignment (L1) is not at least the length of the N2 gene + 2xtrim_spacer, then drop it
# Remove rowIDs of diff dataframe from nucmer_longest_jumpRemoved, and then rbind(rows) of filtered diff column that only contains longest, most syntenic alignment
diff_info <- diff %>%
  dplyr::mutate(n2_gene_len = (end_gene - start_gene)) %>%
  dplyr::mutate(n2_gene_middle = start_gene + (n2_gene_len / 2)) %>%
  dplyr::mutate(slope = ((Et2 - St2) / (end_aln - start_aln))) %>%
  dplyr::mutate(intercept = St2 - (slope * start_aln)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_n2_middleGene = ((slope * n2_gene_middle) + intercept)) %>%
  dplyr::arrange(St2) %>%
  dplyr::group_by(n2_gene, strain) %>% 
  dplyr::mutate(local_dup = ifelse((min(Et2) < max(St2)) & (max(St2) > min(Et2)) & ( (L1 == min(L1) & start_aln < (start_gene - 0) & end_aln > (end_gene + 0)) ), TRUE, FALSE)) %>% 
  dplyr::mutate(WS_n2_middleGene_diff = max(WS_n2_middleGene) - min(WS_n2_middleGene)) %>% # want to include this data, because we probably don't want to keep massive jumps in duplication coordinates
  dplyr::ungroup() 


dist <- ggplot(diff_info) +
  scale_y_continuous(expand = c(0.005, 0)) +
  scale_x_continuous(expand = c(0.005,0)) +
  geom_point(data = diff_info %>% dplyr::filter(local_dup == FALSE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_point(data = diff_info %>% dplyr::filter(local_dup == TRUE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_vline(xintercept = 100, color = "gray30", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'gray', alpha = 0.008) +
  geom_hline(yintercept = 0.05, color = "gray30", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05, fill = 'gray', alpha = 0.008) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.title = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 14, color = 'black', face = 'bold'),
    panel.border = element_rect(fill = NA)) +
  labs(x = "WS alignment coordinate difference at N2 gene locus (kb)", y = "Proportional difference in size of WS alignments (1 - (min(L2) / max(L2)))")
dist

# How many are local (contained in longer WS contig alignment) duplications?
local_dup_counts <- diff_info %>%
  dplyr::count(local_dup) 
# 175 are local, or should be kept due to syntenic context with the threshold of having to span the N2 gene plus 6kb upstream and downstream
# 10,449 are local, or should be kept due to syntenic context when I lower the threshold of synteny to only having to span the N2 gene, and not 6kb up- and downstream

gene_jump_dist <- ggplot(data = diff_info) +
  geom_histogram(aes(x = WS_n2_middleGene_diff / 1e3), bins = 200, fill = 'gray30') +
  scale_y_continuous(expand = c(0.001, 0), labels = scales::label_number(scale = 1e-3)) +
  scale_x_continuous(expand = c(0.001,0)) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color ='black'),
    # plot.margin = margin(r = 10, t = 10, l = 10),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank()) +
  ylab("Count (thousands)")
gene_jump_dist

L2_diff_dist <- ggplot(data = diff_info) +
  geom_histogram(aes( y= diff_fraction), bins = 200, fill = 'gray30') +
  scale_y_continuous(expand = c(0.001, 0), position = "right") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-3)) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = 'black'),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16, color = 'black'),
    # plot.margin = margin(l = 20, t = 10, r = 10, b = 23),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.bottom = element_blank()) +
  xlab("Count (thousands)")
L2_diff_dist


gene_jump_dist_clean <- gene_jump_dist + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

L2_diff_dist_clean <- L2_diff_dist + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

top_row <- plot_grid(gene_jump_dist_clean, NULL, ncol = 2, rel_widths = c(0.8, 0.20))
middle_row <- plot_grid(dist, L2_diff_dist_clean, ncol = 2, rel_widths = c(0.8, 0.20)) #+ theme(plot.margin = margin(t = 20))

final_plot <- plot_grid(top_row, middle_row, nrow = 2, rel_heights = c(0.2, 0.8))

final_plot


# CRAZY EXAMPLE of gene duplication over 10Mb away on the same contig
lkj <- diff_info %>% dplyr::filter(WS_n2_middleGene_diff > 10100000 & diff_fraction > 0.75 & local_dup == T & n2_gene == "WBGene00013102")
lkjplt <- ggplot(diff_info %>% dplyr::filter(n2_gene == "WBGene00013102" & strain == "ECA2187")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.5) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6), color = 'blue', linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    axis.text = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
  labs(y = "ECA2187 contig position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00013102 : IV")
# coord_cartesian(xlim = c(7.4,7.408), ylim = c(3,3.05))
lkjplt


# asdf <- diff_info %>% dplyr::filter(WS_n2_middleGene_diff > 400000 & diff_fraction > 0.75 & local_dup == T)
# boo <- ggplot(diff_info %>% dplyr::filter(n2_gene == "WBGene00012012" & strain == "ECA1887")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA1887 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00012012 : IV")
# # coord_cartesian(xlim = c(7.4,7.408), ylim = c(3,3.05))
# boo
# 
# 
# 
# pltExL2 <- diff_info %>% dplyr::filter(n2_gene == "WBGene00000638") %>% dplyr::filter(strain == "ECA2607")
# dfFilt_gene <- ggplot(pltExL2) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA2607 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00000638 : I")
#   # coord_cartesian(xlim = c(7.4,7.408), ylim = c(3,3.05))
# dfFilt_gene



 
# Example of two alignments of the same length from the same contig that are very far apart
# pltExL3 <- diff_info %>% dplyr::filter(n2_gene == "WBGene00012066") %>% dplyr::filter(strain == "CB4856")
# dfFilt_gene2 <- ggplot(pltExL3) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "CB4856 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00012066 : V")
# dfFilt_gene2

# second largest on x-axis
# pltExL4 <- diff_info %>% dplyr::filter(n2_gene == "WBGene00010874" & strain == "CB4852")
# dfFilt_gene4 <- ggplot(pltExL4) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "CB4852 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00010874 : I")
# dfFilt_gene4
# 
# super small on x-axis but large (ish) on y-axis
# pltExL5 <- diff_info %>% dplyr::filter(strain == "ECA1943") %>% dplyr::filter(n2_gene =="WBGene00011719")
# dfFilt_gene8 <- ggplot(pltExL5) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA1943 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00011719 : IV")
# dfFilt_gene8
# 
# # diff value close to max and larger of the two points ~370kb on x-axis
# pltExL5 <- diff_info %>% dplyr::filter(n2_gene == "WBGene00010874" & strain == "NIC2")
# dfFilt_gene3 <- ggplot(pltExL5) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "NIC2 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00010874 : II") 
# dfFilt_gene3
# 
# # high diff value and right near WS gene coordinate jump cutoff of 100kb
# diff_info_100kb <- diff_info %>% dplyr::filter(WS_n2_middleGene_diff < 100000 & diff > 250000)
# pltExL5 <- diff_info_100kb %>% dplyr::filter(n2_gene == "WBGene00001937" & strain == "ECA738")
# dfFilt_gene4 <- ggplot(pltExL5) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA738 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00001937 : IV") 
# dfFilt_gene4
# 
# # trying to find a diff cutoff
# diff_info_diff_cutoffparsing <- diff_info %>% dplyr::filter(WS_n2_middleGene_diff > 100000 & (diff_fraction < 0.06 & diff_fraction > 0.04)) %>% dplyr::arrange(diff_fraction,WS_n2_middleGene_diff)
# # WBGene00005461 - ECA2377 : keep  
# # WBGene00020652 - ECA742 : keep 
# # WBGene00000875 - LCK34 : keep - 188kb
# # WBGene00001686 - ECA742 : DON'T KEEP - 199kb
# # WBGene00010874 - ECA706 : DON'T KEEP - 100kb
# 
# # >20 is where is really starts to fall off and the alignments are barely expanding past the N2 gene coordinates
# 
# pltExL6 <- diff_info_diff_cutoffparsing %>% dplyr::filter(n2_gene == "WBGene00003466" & strain == "CB4856")
# dfFilt_gene5 <- ggplot(pltExL6) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "CB4856 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00003466 : IV")
# dfFilt_gene5

# 
# 
# plotdf <- diff_info %>% dplyr::filter(n2_gene == "WBGene00022590") %>% dplyr::filter(strain == "AB1")
# dfFilt_gene <- ggplot(plotdf) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   geom_vline(xintercept = plotdf$n2_gene_middle / 1e6, color = 'blue') +
#   geom_hline(yintercept = plotdf$WS_n2_middleGene / 1e6, color = 'blue') +
#   geom_text(data = plotdf %>% dplyr::filter(L2 == max(L2)), aes(x = ((n2_gene_middle / 1e6) - 0.07), y = ((WS_n2_middleGene / 1e6) + 0.03), label = paste0("Difference in WS coordinates: ", round(WS_n2_middleGene_diff))), vjust = -1, color = 'blue', fontface = 'bold', size = 5) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "AB1 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00022590 : V")
# dfFilt_gene


diff_filtered <- diff_info %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(
    drop_smaller = if (
      any(WS_n2_middleGene_diff > 100000 & diff_fraction > 0.05 & local_dup == FALSE)
    ) {
      (L1 != max(L1, na.rm = TRUE)) & (local_dup == FALSE)
    } else {
      FALSE  # keep both
    }
  ) %>%
  dplyr::filter(!drop_smaller) %>%
  dplyr::ungroup()

diff_filtered_twos <- diff_filtered %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(dplyr::n() == 2) %>%
  dplyr::ungroup()
  
# dist_filt <- ggplot(diff_filtered_twos) +
#   scale_y_continuous(expand = c(0.005, 0)) +
#   scale_x_continuous(expand = c(0.005,0)) +
#   geom_point(data = diff_filtered_twos %>% dplyr::filter(local_dup == FALSE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
#   geom_point(data = diff_filtered_twos %>% dplyr::filter(local_dup == TRUE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
#   geom_vline(xintercept = 100, color = "gray30", size = 2, linetype="dashed") +
#   geom_rect(xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'gray', alpha = 0.008) +
#   geom_hline(yintercept = 0.05, color = "gray30", size = 2, linetype="dashed") +
#   geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05, fill = 'gray', alpha = 0.008) +
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     # axis.title = element_blank(),
#     axis.text = element_text(size = 14, color = 'black'),
#     axis.title = element_text(size = 14, color = 'black', face = 'bold'),
#     panel.border = element_rect(fill = NA)) +
#   labs(x = "WS alignment coordinate difference at N2 gene locus (kb)", y = "Proportional difference in size of WS alignments (1 - (min(L2) / max(L2)))")
# dist_filt


# testplt <- diff_filtered %>% dplyr::filter(n2_gene == "WBGene00010874" & strain == "CB4852") 
# dfFilt_gene4 <- ggplot(testplt) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "CB4852 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00010874 : I")
# dfFilt_gene4

# Remove n2gene-strain pairs with two alignments from nucmer_longest_jumpRemoved, and then append filtered alignments
rows_to_remove <- diff$row_id

nucmer_longest_jumpRemoved_doublesGone <- nucmer_longest_jumpRemoved %>%
  dplyr::filter(!row_id %in% rows_to_remove)

add_to_nucmerLongest <- diff_filtered %>%
  dplyr::select(chrom,start_aln,end_aln,L1,longest_contig,WSS,WSE,L2,LENQ,strain,start_gene,end_gene,n2_gene,N2,nalign,ntig,tigsize,inv,St2,Et2,leadDiff,jump,run_id,len)

nucmer_longest_jumpRemoved_updated <- nucmer_longest_jumpRemoved_doublesGone %>%
  dplyr::bind_rows(add_to_nucmerLongest)


# okay <- diff_filtered %>%
#   dplyr::filter(local_dup == TRUE & diff_fraction > 0.05 & WS_n2_middleGene_diff > 100000)
# # A smaller alignment that we would want to keep
# yup <- diff_filtered %>% dplyr::filter(n2_gene == "WBGene00021707") %>% dplyr::filter(strain == "ECA923")  #%>% dplyr::filter(L1 == min(L1))
# yupidk <- ggplot(yup) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA923 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00021707 : V")
# # coord_cartesian(xlim = c(7.4,7.408), ylim = c(3,3.05))
# yupidk

# ======================================================================================================================================================================================== #
# Trimming alignment(s) to ROI (N2 gene) 
# ======================================================================================================================================================================================== #
# # ============= use trim spacer on dataset needed for faceting extreme CNVs genes and plotting trimmed alignments ========================== #
subset <- nucmer_longest_jumpRemoved_updated %>%
  dplyr::filter(n2_gene == "WBGene00044801" | n2_gene == "WBGene00011254" | n2_gene == "WBGene00022486")

trim_spacer = 5e3 # trimming to 5kb on either side of the N2 gene
tigTrim_subset <- subset %>%
  dplyr::arrange(n2_gene,strain,start_aln) %>%
  dplyr::mutate(unchanged_start_aln = start_aln, unchanged_end_aln = end_aln) %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scale_distortion = ((L2 - L1)/L1), # get the distortion of the WS to N2 coordinate transformation - is the slope one?
    rboundDist = max(start_aln, end_aln) - end_gene,
    WSE = ifelse(rboundDist > trim_spacer & inv == F, WSE + (round(scale_distortion*rboundDist)) - (rboundDist - trim_spacer), WSE),
    WSE = ifelse(rboundDist > trim_spacer & inv == T, WSE - (round(scale_distortion*rboundDist)) + (rboundDist - trim_spacer), WSE),
    end_aln = ifelse(rboundDist > trim_spacer,end_aln - (rboundDist - trim_spacer), end_aln),
    lboundDist = start_gene - min(start_aln, end_aln),
    WSS = ifelse(lboundDist > trim_spacer & inv == F, WSS + (round(scale_distortion*lboundDist)) + (lboundDist - trim_spacer), WSS),
    WSS = ifelse(lboundDist > trim_spacer & inv == T, WSS - (round(scale_distortion*lboundDist)) - (lboundDist - trim_spacer), WSS),
    start_aln = ifelse(lboundDist > trim_spacer, start_aln + (lboundDist - trim_spacer), start_aln)) %>%
  dplyr::ungroup()

gene1 <- tigTrim_subset %>% dplyr::filter(n2_gene == "WBGene00011254")
g7plt <- ggplot(gene1) +
  geom_rect(data = gene1 %>% dplyr::distinct(strain, .keep_all = T), aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.5) +
  facet_wrap(~strain, scales = "free") + #strip.position = "right") +
  # geom_segment(data = nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00011254"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'black', linewidth = 1) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = "blue", linewidth = 1) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle(expression(italic("R12G8.1") * ": V"))
g7plt


examination <- tigTrim_subset %>% dplyr::filter(n2_gene == "WBGene00022486" & strain == "JU1581")
explot <- ggplot(examination) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.5) +
  # facet_wrap(~strain, scales = "free") + #strip.position = "right") +
  # geom_segment(data = nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00011254"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'black', linewidth = 1) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = "blue", linewidth = 1) +
  theme_bw() 
explot



# ========================================================== #



trim_spacer = 5e3 # trimming to 5kb on either side of the N2 gene
tigTrim <- nucmer_longest_jumpRemoved_updated %>%
  dplyr::arrange(n2_gene,strain,start_aln) %>% 
  dplyr::mutate(unchanged_start_aln = start_aln, unchanged_end_aln = end_aln) %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(
    scale_distortion = ((L2 - L1)/L1), # get the distortion of the WS to N2 coordinate transformation - is the slope one?
    rboundDist = max(start_aln, end_aln) - end_gene,
    Et2 = ifelse(rboundDist > trim_spacer, Et2 + (round(scale_distortion*rboundDist)) - (rboundDist - trim_spacer), Et2),
    end_aln = ifelse(rboundDist > trim_spacer,end_aln - (rboundDist - trim_spacer),end_aln),
    lboundDist = start_gene - min(start_aln, end_aln),
    St2 = ifelse(lboundDist > trim_spacer, St2 + (round(scale_distortion*lboundDist)) + (lboundDist - trim_spacer), St2),
    start_aln = ifelse(lboundDist > trim_spacer,start_aln + (lboundDist - trim_spacer), start_aln)) %>%
  dplyr::ungroup()
 

# Confirming tigTrim is working correctly and not distorting coordinates
# g10plt <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00044801")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   facet_wrap(~strain, scales = "free", strip.position = "right") +
#   geom_segment(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00044801"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6), color = 'black') +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig)) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00044801")
# g10plt
# 
# # ^ Same as plot above but for WSS and WSE
# g11plt <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00044801")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   facet_wrap(~strain, scales = "free", strip.position = "right") +
#   geom_segment(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00044801"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'black', size = 3) +
#   geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), size = 1) +
#   theme_bw() +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00044801")
# g11plt


# # Example of positive scale_distortion
# ctg_trim <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00013300") %>% dplyr::filter(strain == "ECA369")) +
#     geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#     geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#     theme(
#       legend.position = 'none',
#       panel.background = element_blank(),
#       panel.grid = element_blank(),
#       panel.border = element_rect(fill = NA),
#       plot.title = element_text(size = 14, color = 'black')) +
#     labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# ctg_trim
# 
# testYUP <- tigTrim %>% dplyr::filter(n2_gene == "WBGene00010874" & strain == "CB4852") 
# trimTig1 <- ggplot(testYUP) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "CB4852 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00010874 : I")
# trimTig1
# 
# 
# ctg_trim4 <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00000638") %>% dplyr::filter(strain == "ECA2607")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA2607 coord", x = "N2 coord", title = "WBGene00000638 : I")
# ctg_trim4



# ======================================================================================================================================================================================== #
# Adding WS genes to each alignment, and then collapsing by N2 gene to get a matrix of syntenic genes #
# ======================================================================================================================================================================================== #
wsg <- data.table::as.data.table(ws_genes)

clean_tigTrim <- tigTrim %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::select(chrom, unchanged_start_aln, unchanged_end_aln, start_aln, end_aln, L1, start_gene, end_gene, n2_gene, WSS, WSE, inv, St2, Et2, longest_contig, nalign, L2, strain) 

filt_nucm_long <- data.table::as.data.table(clean_tigTrim)

data.table::setnames(filt_nucm_long, c("longest_contig", "St2", "Et2"), c("contig", "start", "end")) # use St2 Et2 because the gff coordinates will be reported in this way

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
  
# write.table(syntelog_matrix, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/syntelog_matrix.tsv", sep = '\t', col.names = T, row.names = F, quote = F)
  
  

# # Confirm this works by plotting N2 gene coordinates in x system and WS gene coordinates in y system
# syn1 <- ggplot(joined %>% dplyr::filter(strain == "AB1", n2_gene == "WBGene00015153")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00015153") %>% dplyr::filter(strain == "AB1"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00015153 : V")
# syn1

# rando <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00012259")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   facet_wrap(~strain, scales = 'free') +
#   # geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00015153") %>% dplyr::filter(strain == "AB1"),
#             # aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00015153 : V")
# rando

# 
# ugh <- ggplot(joined %>% dplyr::filter(strain == "ECA1208", n2_gene == "WBGene00011875")) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.15) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.25) +
#   geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6),color = 'blue', linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00011875" & strain == "ECA1208"),
#             aes(x = 10.212, y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, size = 5, color = '#3B2F2F', fontface = 'bold') +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00011875" & strain == "ECA1208" & WSS == "4936031"), 
#             aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = 4.926, label = "N2 gene: WBGene00020062"), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
#   
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     axis.title = element_text(size = 14, color = 'black', face = 'bold'),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 14, color = 'black'),
#     panel.border = element_rect(fill = NA)) +
#     # plot.title = element_text(size = 14, color = 'black', hjust = 0.5, face = 'bold')) +
#   labs(y = "ECA1208 coord", x = "N2 coord")
# ugh
# 
# # There is only one alignment that has a predicted gene - resolves itself
# ugh2 <- ggplot(joined %>% dplyr::filter(strain == "ECA768", n2_gene == "WBGene00000864")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00000864") %>% dplyr::filter(strain == "ECA768"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA768 coord", x = "N2 coord", title = "WBGene00000864 : II")
# ugh2
 
# hmm <- ggplot(clean_tigTrim %>% dplyr::filter(n2_gene == "WBGene00000864") %>% dplyr::filter(strain == "ECA768")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 14, color = 'black'),
#     plot.title = element_text(size = 18, hjust = 0.5, color = 'black', face = 'bold')) +
#   labs(y = "ECA768 genome position (Mb)", x = "N2 genome position (Mb)", title = "WBGene00000864 : V")
# hmm
# 
# ugh3 <- ggplot(joined %>% dplyr::filter(strain == "CB4856", n2_gene == "WBGene00012066")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00012066") %>% dplyr::filter(strain == "CB4856"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "CB4856 coord", x = "N2 coord", title = "WBGene00012066 : V")
# ugh3
# 
# 
# synner <- joined %>% dplyr::filter(strain == "ECA923", n2_gene == "WBGene00021707")
# syn0 <- ggplot(synner) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.7) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6), color = "blue", linewidth = 1) +
#   # geom_text(data = syn,
#             # aes(x =  (start_aln + 400) / 1e6 , y = (i.end  - (i.end - i.start + 225) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
#   # geom_text(data = syn,
#             # aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = WSS / 1e6, label = "N2 gene: WBGene00021707"), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     axis.title = element_text(size = 14, color = 'black', face = 'bold'),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 14, color = 'black'),
#     panel.border = element_rect(fill = NA),
#   plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
#   labs(y = "ECA923 genome coordinates (Mb)", x = "N2 genome coordinates (Mb)", title = "WBGene00021707 : V")
# syn0
# 
# 
# syn2 <- ggplot(joined %>% dplyr::filter(strain == "ECA369", n2_gene == "WBGene00013300")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00013300") %>% dplyr::filter(strain == "ECA369"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# syn2
# 
# 
# syn3 <- ggplot(joined %>% dplyr::filter(strain == "ECA2529", n2_gene == "WBGene00020754")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00020754") %>% dplyr::filter(strain == "ECA2529"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA2529 coord", x = "N2 coord", title = "WBGene00020754 : V")
# syn3
# 
# 
# syn4 <- ggplot(joined %>% dplyr::filter(strain == "ECA1202", n2_gene == "WBGene00044801")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00044801") %>% dplyr::filter(strain == "ECA1202"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "ECA1202 coord", x = "N2 coord", title = "WBGene00044801 : V")
# syn4
# 
# syn5 <- ggplot(joined %>% dplyr::filter(strain == "CB4856", n2_gene == "WBGene00009284")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00009284") %>% dplyr::filter(strain == "CB4856"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "CB4856 coord", x = "N2 coord", title = "WBGene00009284 : V")
# syn5
# 
# syntest <- ggplot(joined %>% dplyr::filter(strain == "AB1", n2_gene == "WBGene00022590")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
#   geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00022590") %>% dplyr::filter(strain == "AB1"),
#             aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00022590 : V")
# syntest

# 
syn <- joined %>% dplyr::filter(strain == "AB1", n2_gene == "WBGene00020062")
syn6 <- ggplot(syn) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.7) +
  geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = "blue", linewidth = 1) +
  geom_text(data = syn,
            aes(x =  (start_aln + 475) / 1e6 , y = (i.end  - (i.end - i.start + 550) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 9, fontface = 'bold') +
  geom_text(data = syn,
            aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = WSS / 1e6, label = "N2 gene: nhr-270"), vjust = -1, color = '#3B2F2F', size = 9, fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_text(size = 20, color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(fill = NA)) +
    # plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
  labs(y = "AB1 contig coordinates (Mb)", x = "N2 genome coordinates (Mb)")
syn6
# 
# syn2 <- joined %>% dplyr::filter(strain == "ECA1887", n2_gene == "WBGene00011257")
# syn123 <- ggplot(syn2) +
#   # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.7) +
#   geom_segment(data = diff_filtered %>% dplyr::filter(strain == "ECA1887" & n2_gene == "WBGene00011257"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6), color = "blue", linewidth = 1) +
#   # geom_text(data = syn2,
#             # aes(x =  (start_aln + 400) / 1e6 , y = (i.end  - (i.end - i.start + 225) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
#   # geom_text(data = syn2,
#             # aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = start / 1e6, label = "N2 gene: npax-3"), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     axis.title = element_text(size = 20, color = 'black', face = 'bold'),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 14, color = 'black'),
#     panel.border = element_rect(fill = NA)) +
#   # plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
#   labs(y = "ECA1887 genome coordinates (Mb)", x = "N2 genome coordinates (Mb)")
# syn123




# # ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/syntenicGene.png", dpi = 600, syn6, width = 14, height = 12)
# 
syn2 <- joined %>% dplyr::filter(strain == "ECA1260", n2_gene == "WBGene00011254")
syn7 <- ggplot(syn2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.15) +
  geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'blue', linewidth = 1) +
  geom_text(data = syn2,
            aes(x =  17.7455, y = ((i.start / 1e6) + (i.end / 1e6)) / 2.0087, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 7, fontface = 'bold') +
  geom_text(data = syn2 %>% dplyr::filter(L2 == "5151"),
            aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = 0.31, label = "N2 gene: R12G8.1"), vjust = -1, color = '#3B2F2F', size = 7, fontface = 'bold') +

  # geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00011254") %>% dplyr::filter(strain == "ECA1260"),
            # aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_text(size = 20, color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(fill = NA)) +
    # plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
  labs(y = "ECA1260 contig coordinates (Mb)", x = "N2 genome coordinates (Mb)")
syn7


bwc <- joined %>% dplyr::filter(strain == "JU1581", n2_gene == "WBGene00022486")
bwc1 <- ggplot(bwc) +
  geom_rect(data = bwc %>% dplyr::filter(attributes == "g22825"), aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.15) +
  geom_vline(xintercept = 1239412 / 1e6, color = 'black', linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = 1241124 / 1e6, color = 'black', linetype = "dashed", size = 1.5) +
  # geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'blue', linewidth = 1) +
  geom_segment(data = nucmer_longest %>% dplyr::filter(strain == "JU1581" & n2_gene == "WBGene00022486"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = "blue", linewidth = 1) +
  geom_text(data = bwc %>% dplyr::filter(attributes == "g22825"),
            aes(x =  1.230, y = (i.start / 1e6) - 0.001, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 7, fontface = 'bold') +
  geom_text(data = bwc, # gunna need to fix dis shit!!
            aes(x =  (end_gene) / 1e6 + 0.001, y = 0.05, label = "N2 gene: fbxa-61"), vjust = -1, color = '#3B2F2F', size = 7, fontface = 'bold') +
  
  # geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00011254") %>% dplyr::filter(strain == "ECA1260"),
  # aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_text(size = 20, color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  # plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
  labs(y = "JU1581 contig coordinates (Mb)", x = "N2 genome coordinates (Mb)")
bwc1


# 
# bwc_yup <- joined %>% dplyr::filter(strain == "ECA1825", n2_gene == "WBGene00044801")
# bwc2 <- ggplot(bwc_yup) +
#   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.15) +
#   geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'blue', linewidth = 1) +
#   # geom_text(data = syn2,
#             # aes(x =  17.7455, y = ((i.start / 1e6) + (i.end / 1e6)) / 2.0087, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 7, fontface = 'bold') +
#   # geom_text(data = syn2 %>% dplyr::filter(L2 == "5151"), # gunna need to fix dis shit!!
#             # aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = 0.31, label = "N2 gene: clc-30"), vjust = -1, color = '#3B2F2F', size = 7, fontface = 'bold') +
#   
#   # geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00011254") %>% dplyr::filter(strain == "ECA1260"),
#   # aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     axis.title = element_text(size = 20, color = 'black', face = 'bold'),
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 14, color = 'black'),
#     panel.border = element_rect(fill = NA)) +
#   # plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
#   labs(y = "ECA1825 contig coordinates (Mb)", x = "N2 genome coordinates (Mb)")
# bwc2



### Example of how plotting the inverted version messes up how the WS gene models are called
oop <- joined %>% dplyr::filter(n2_gene == "WBGene00011254"& strain == "JU2617")
oops <- ggplot(oop) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.end / 1e6, ymax = i.start / 1e6), fill = "#DB6333", alpha = 0.3) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.15) +
  geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'blue', linewidth = 1) +
  # The non-inverted, trimmed alignment
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6), color = 'blue', linewidth = 1) +
  facet_wrap(~strain, scales = "free") +
  # geom_text(data = syn2,
            # aes(x =  17.745, y = ((i.start / 1e6) + (i.end / 1e6)) / 2.0083, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
  # geom_text(data = syn2 %>% dplyr::filter(L2 == "5151"),
            # aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = 0.31, label = "N2 gene: WBGene00011254"), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +

  # geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00011254") %>% dplyr::filter(strain == "ECA1260"),
  # aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_blank(),
    # axis.ticks = element_blank(),
    panel.grid = element_blank(),
    # axis.text = element_blank(),
    panel.border = element_rect(fill = NA))
  # plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'bold')) +
  # labs(y = "ECA1260 genome coordinates (Mb)", x = "N2 genome coordinates (Mb)")
oops










# output_list <- list()
# for i in 1:length(trimmed_final){
#   group_list <- list()
#   for j in 1:nrow(trimmed_final[[i]]){
#     dplyr::filter(ws_genes$strain == strain & ws_genes$contig == contig) %>%
#       dplyr::filter(ws_genes$start >= start & ws_genes$end <= end)
#
#     group_list[[j]] <- list(n2_gene, ws_genes$contig, ws_genes$start, ws_genes$end, ws_genes$attributes)
#   }
#   output_list[[i]] <- ldply(group_list[[j]],data.frame) %>% pivot_wider(names_from = ws_genes$attributes, values_from = syntenic_gene)
# }