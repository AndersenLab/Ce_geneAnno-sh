library(dplyr)
library(tidyr)
library(ggplot2)
library(fuzzyjoin)
library(purrr)
library(stringr)
library(data.table)

# ======================================================================================================================================================================================== #
# Loading N0.tsv with transcripts converted to genes and plotting gene set classification #
# ======================================================================================================================================================================================== #
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0504_genes.tsv")

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

all_relations <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))


#### Plotting ####
private_freq = (1/(length(strainCol_c2_u)))

classification <- all_relations %>%
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
  dplyr::count(freq, class) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) 


gs_allOrtho <- ggplot(data = classification, aes(x = freq * 100, y = percent, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("HOGs") + # longest isoform
  xlab("Frequency") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=16, color = 'black'),
    axis.text = element_text(size=14, color = 'black')
  )
gs_allOrtho

# HOG_class_count <- classification %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarise(n_HOG = sum(n)) %>%
#   dplyr::ungroup()



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
ortho_genes_dd_rowid <- ortho_genes_dd %>% dplyr::rename(rowid = HOG)

# Extract rowids for simple and complex HOGs
simple_rowids <- all_relations_rowid %>% 
  dplyr::filter(if_all(2:ncol(.), ~ is.na(.) | . <= 1)) %>% 
  dplyr::pull(rowid)

complex_rowids <- all_relations_rowid %>% 
  dplyr::filter(if_any(2:ncol(.), ~ . > 1)) %>% 
  dplyr::pull(rowid)

simple_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% simple_rowids)
complex_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% complex_rowids) ### further decompress into only complex orthogroups that have ONE N2 gene

n2_goi <- complex_HOGS%>%
  dplyr::select(N2) %>%
  tidyr::separate_rows(N2, sep = ",\\s*") %>%
  dplyr::filter(!is.na(N2)) 

n2_genes <- all_genes_strain %>% dplyr::filter(strain == "N2") %>% dplyr::rename(chrom = contig) %>%
  dplyr::filter(reduce(n2_goi$N2, ~ .x | str_detect(attributes, fixed(.y)), .init = FALSE))

# ======================================================================================================================================================================================== #
# Extracting the longest WS contig alignment for every N2 gene coordinate #
# ======================================================================================================================================================================================== #
nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, contig, WSS, WSE, L2, LENQ, strain)

# Join where an N2 gene is fully contained within an alignment
n2_genes_align <- fuzzyjoin::interval_inner_join(
  x = n2_genes,
  y = nucmer_ranges,
  by = c("start", "end"),
  type = "within" 
) %>%
  dplyr::filter(chrom.x == chrom.y) %>%
  dplyr::rename(start_gene = start.x, end_gene = end.x, start_aln = start.y, end_aln = end.y, n2_gene = attributes, strain = strain.y, N2 = strain.x)


nucmer_longest <- n2_genes_align %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, contig, strain) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L2)) %>% #summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times... does this mean it is a "better" alignment?
  dplyr::ungroup() %>% ##### ADD AN ARGUMENT FOR IF TWO CONTIGS HAVE THE SAME ALINGMENT, TO CHOOSE THE LONGEST LENQ
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()

# #### Testing how selection of "best" contig works with dot plots ####
# all_ctg_align <- n2_genes_align %>%
#   dplyr::group_by(strain, n2_gene) %>%
#   dplyr::mutate(nalign = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(n2_gene,contig, strain) %>%
#   dplyr::mutate(ntig= n()) %>%
#   dplyr::mutate(tigsize=sum(L2)) %>%
#   dplyr::ungroup() 

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
# g7plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00011254")) +
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
#         plot.title = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5)) +
#   ggtitle("WBGene00011254")
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
# g10plt <- ggplot(onekb %>% dplyr::filter(n2_gene == "WBGene00044801")) +
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
# g10plt


# ======================================================================================================================================================================================== #
# Removing large jumps in alignment #
# ======================================================================================================================================================================================== #transformed_coords <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/nucmer_runs/115_WI_transformed_coords_FIXED.tsv", col_names = F) 
nucmer_longest_jump <- nucmer_longest %>%
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) %>% # addresses inverted alignments for tigTrim
  dplyr::arrange(St2) %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(leadDiff = lead(St2)-Et2)
                
# ggplot(nucmer_longest_jump) +
#   geom_histogram(aes(x = leadDiff), bins = 300) +
#   geom_vline(xintercept = 4.5E5) +
#   geom_vline(xintercept = -4.5E5) +
#   # geom_vline(xintercept = -1.5E5) +
#   # coord_cartesian(x = c(-300000, 300000)) + 
#   #scale_y_log10() +
#   theme_bw() 
#   # xlim(-6E5,6E5)


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
  
  




# ======================================================================================================================================================================================== #
# Selecting longest alignment for contigs that have multiple alignments of the same length after filtering for jumps 
# & trimming alignments to ROI (N2 gene) to ensure synteny 
# ======================================================================================================================================================================================== #
# Extracting situations where there are two alignments of a contig, neither are removed from the jump, and assessing their pairwise difference in alignment length
diff <- nucmer_longest_jumpRemoved %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(dplyr::n() == 2) %>%
  dplyr::arrange(L2, .by_group = TRUE) %>%
  dplyr::mutate(diff = abs(diff(L2)[1])) %>% 
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
  # geom_vline(xintercept = 0.5, color = 'blue', size = 2) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  xlab("Difference in alignment jumps (absolute value) (Mb)")
diff_dist_plotjumps      
              
negJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00009423") %>% dplyr::filter(strain == "JU310")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "JU310 coord", x = "N2 coord", title = "WBGene00009423 : V")
negJump


zeroJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00000680") %>% dplyr::filter(strain == "PX179")) +
    geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
    geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA),
      plot.title = element_text(size = 14, color = 'black')) +
    labs(y = "PX179 coord", x = "N2 coord", title = "WBGene00000680 : I")
zeroJump

lrgJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00000638") %>% dplyr::filter(strain == "XZ1513")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "XZ1513 coord", x = "N2 coord", title = "WBGene00000638 : I")
lrgJump

# TWO DIFFERENT ALIGNMENTS FOR THE SAME CONTIG THAT ARE THE SAME SIZE
largestJump <- ggplot(diff %>% dplyr::filter(n2_gene == "WBGene00012066") %>% dplyr::filter(strain == "CB4856")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "CB4856 coord", x = "N2 coord", title = "WBGene00018439 : V")
largestJump
  
contigDup <- ggplot() +
  geom_segment(data = nucmer %>% dplyr::filter(contig == "ptg000004l") %>% dplyr::filter(strain == "CB4856"), aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = diff %>% dplyr::filter(strain == "CB4856" & n2_gene == "WBGene00012066"), aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "green") +
  coord_cartesian(xlim = c(19.2,19.45)) + 
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "CB4856 coord", x = "N2 coord", title = "CB4856 ptg000004l")
contigDup

# Removing alignments that are VERY small and distant from "main" alignment
# Need to find the WS coordinates that overlap with the N2 gene for each alignment, then calculate the mean (middle) WS coordinate, 
# and then if the difference is >trim_spacer (need to account for scale_distortion between WS and N2) and the smaller alignment (L1) is not at least the length of the N2 gene + 2xtrim_spacer, then drop it
# Remove rowIDs of diff dataframe from nucmer_longest_jumpRemoved, and then rbind(rows) of filtered diff column that only contains longest, most syntenic alignment
trim_spacer = 5e3 # used later in tigTrim

diff_filtered <- diff %>%
  dplyr::mutate(n2_gene_len = (end_gene - start_gene)) %>%
  dplyr::mutate(n2_gene_middle = start_gene + (n2_gene_len / 2)) %>%
  dplyr::mutate(slope = ((Et2 - St2) / (end_aln - start_aln))) %>%
  dplyr::mutate(intercept = St2 - (slope * start_aln)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_n2_middleGene = ((slope * n2_gene_middle) + intercept)) %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(WS_n2_middleGene_diff = max(WS_n2_middleGene) - min(WS_n2_middleGene)) %>%
  dplyr::rowwise()
  dplyr::mutate(distAlignGreater = ifelse((lead(St2) > Et2 & lead(Et2) > Et2) | (lead(St2) & lead(Et2) < St2), TRUE, FALSE)) %>% #NEED TO HAVE LEAD(ST2) AND LEAD(E2) BE GREATER Need both coordinates to be larger???
  dplyr::ungroup() %>%
  dplyr::mutate(WS_jump_greaterTrimSpacer = ifelse((abs((((L2 - L1)/L1)) * WS_n2_middleGene_diff) + WS_n2_middleGene_diff) > trim_spacer, TRUE, FALSE)) %>% # account for WS coordinate transformation and add/remove to reflect length of trim_spacer in N2 coordinates 
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(L2 == ifelse((WS_jump_greaterTrimSpacer == T & (L1 < (n2_gene_len + 10000))), max(L2),L2)) %>% # filtering if the jump in WS coordinates is greater than trim_spacer, and the L1 alignment is less than n2_gene + trim_spacer on both sides
  dplyr::ungroup()

# dplyr::select(chrom.x,type,start_gene,end_gene,strand,n2_gene,N2,chrom.y,start_aln,end_aln,L1,
#               longest_contig,WSS,WSE,L2,strain,nalign,ntig,tigsize,inv,St2,Et2,leadDiff,jump,run_id,len,sumlen) # only selecting columns that are in nucmer_longest_jumpRemoved for rbinding back together

plotdf <- diff_filtered %>% dplyr::filter(n2_gene == "WBGene00022590") %>% dplyr::filter(strain == "AB1")
dfFilt_gene <- ggplot(plotdf) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  geom_vline(xintercept = plotdf$n2_gene_middle / 1e6, color = 'blue') +
  geom_hline(yintercept = plotdf$WS_n2_middleGene / 1e6, color = 'blue') +
  geom_text(data = plotdf %>% dplyr::filter(L2 == max(L2)), aes(x = ((n2_gene_middle / 1e6) - 0.07), y = ((WS_n2_middleGene / 1e6) + 0.03), label = paste0("Diff in WS coordinates: ", round(WS_n2_middleGene_diff))), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00022590 : V")
dfFilt_gene




# lngctg <- ggplot(lng_ctg %>% dplyr::filter(n2_gene == "WBGene00013300") %>% dplyr::filter(strain == "ECA369")) +
#     geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#     geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#     theme(
#       legend.position = 'none',
#       panel.background = element_blank(),
#       panel.grid = element_blank(),
#       panel.border = element_rect(fill = NA),
#       plot.title = element_text(size = 14, color = 'black')) +
#     labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# lngctg
# 
# lngctg2 <- ggplot(lng_ctg %>% dplyr::filter(n2_gene == "WBGene00044801") %>% dplyr::filter(strain == "ECA1202")) +
#   geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
#   geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 14, color = 'black')) +
#   labs(y = "WBGene00044801 coord", x = "N2 coord", title = "WBGene00013300 : IV")
# lngctg2


trim_spacer = 5e3 # trimming to 5kb on either side
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
tigTrim <- nucmer_longest_jumpRemoved %>%
  dplyr::arrange(n2_gene,strain,start_aln) %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scale_distortion = ((L2 - L1)/L1), # get the distortion of the WS to N2 coordinate transformation - is the slope one?
    rboundDist = max(start_aln, end_aln) - end_gene,
    Et2 = ifelse(rboundDist > trim_spacer,Et2 + (scale_distortion*rboundDist) - (rboundDist - trim_spacer),Et2), # correctly remove 
    end_aln = ifelse(rboundDist > trim_spacer,end_aln - (rboundDist - trim_spacer),end_aln),
    lboundDist = start_gene - min(start_aln, end_aln),
    St2 = ifelse(lboundDist > trim_spacer,St2 + (scale_distortion*lboundDist) + (lboundDist - trim_spacer),St2),
    start_aln = ifelse(lboundDist > trim_spacer,start_aln + (lboundDist - trim_spacer),start_aln))


# oneex <- nucmer_longest %>% dplyr::filter(n2_gene == "WBGene00013300" & strain == "ECA369")
# twoex <- nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00013300" & strain == "ECA369")
# threeex <- tigTrim %>% dplyr::filter(n2_gene == "WBGene00013300" & strain == "ECA369")

 
# Example of positive scale_distortion
ctg_trim <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00013300") %>% dplyr::filter(strain == "ECA369")) +
    geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
    geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
    theme(
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA),
      plot.title = element_text(size = 14, color = 'black')) +
    labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
ctg_trim

# Example of negative scale_distortion
ctg_trim1 <- ggplot(nucmer_longest_jumpRemoved %>% dplyr::filter(n2_gene == "WBGene00269422") %>% dplyr::filter(strain == "ECA2111")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "JU311 coord", x = "N2 coord", title = "WBGene00269422 : X")
ctg_trim1

ctg_trim2 <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00269422") %>% dplyr::filter(strain == "ECA2111")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "ECA2111 coord", x = "N2 coord", title = "WBGene00269422 : X")
ctg_trim2


ctg_trim3 <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00022590") %>% dplyr::filter(strain == "AB1")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00022590 : V")
ctg_trim3

ctg_trim4 <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00000638") %>% dplyr::filter(strain == "ECA2607")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = St2 / 1e6, yend = Et2 / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) +
  labs(y = "ECA2607 coord", x = "N2 coord", title = "WBGene00000638 : I")
ctg_trim4




# ======================================================================================================================================================================================== #
# Adding WS genes to each alignment, and then collapsing by N2 gene to get a matrix of syntenic genes #
# ======================================================================================================================================================================================== #
wsg <- data.table::as.data.table(ws_genes)

clean_tigTrim <- tigTrim %>%
  dplyr::select(chrom.x, start_aln, end_aln, L1, sstart_gene, end_gene, n2_gene, WSS, WSE, inv, St2, Et2, longest_contig, L2, strain) %>%
  dplyr::rename(chrom = chrom.x)

filt_nucm_long <- data.table::as.data.table(clean_tigTrim)

data.table::setnames(filt_nucm_long, c("longest_contig", "St2", "Et2"), c("contig", "start", "end")) #use St2 Et2 because the gff coordinates will be reported in this way

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
  )
  

# Confirm this works by plotting N2 gene coordinates in x system and WS gene coordinates in y system
syn1 <- ggplot(joined %>% dplyr::filter(strain == "AB1", n2_gene == "WBGene00015153")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00015153") %>% dplyr::filter(strain == "AB1"), 
            aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00015153 : V")
syn1

syn2 <- ggplot(joined %>% dplyr::filter(strain == "ECA369", n2_gene == "WBGene00013300")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00013300") %>% dplyr::filter(strain == "ECA369"), 
            aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "ECA369 coord", x = "N2 coord", title = "WBGene00013300 : IV")
syn2


syn3 <- ggplot(joined %>% dplyr::filter(strain == "ECA2529", n2_gene == "WBGene00020754")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00020754") %>% dplyr::filter(strain == "ECA2529"), 
            aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "ECA2529 coord", x = "N2 coord", title = "WBGene00020754 : V")
syn3


syn4 <- ggplot(joined %>% dplyr::filter(strain == "ECA1202", n2_gene == "WBGene00044801")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00044801") %>% dplyr::filter(strain == "ECA1202"), 
            aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "ECA1202 coord", x = "N2 coord", title = "WBGene00044801 : V")
syn4

syn5 <- ggplot(joined %>% dplyr::filter(strain == "CB4856", n2_gene == "WBGene00009284")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00009284") %>% dplyr::filter(strain == "CB4856"), 
            aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "CB4856 coord", x = "N2 coord", title = "WBGene00009284 : V")
syn5

syntest <- ggplot(joined %>% dplyr::filter(strain == "AB1", n2_gene == "WBGene00022590")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  geom_text(data = joined %>% dplyr::filter(n2_gene == "WBGene00022590") %>% dplyr::filter(strain == "AB1"), 
            aes(x =  start_gene / 1e6 , y = (i.end  - (i.end - i.start) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = 'blue', fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "AB1 coord", x = "N2 coord", title = "WBGene00022590 : V")
syntest




syn <- joined %>% dplyr::filter(strain == "AB1", n2_gene == "WBGene00020062")
syn6 <- ggplot(syn) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.7) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6), color = "blue", linewidth = 1) +
  geom_text(data = syn, 
            aes(x =  (start_aln + 400) / 1e6 , y = (i.end  - (i.end - i.start + 225) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
  geom_text(data = syn, 
            aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = WSS / 1e6, label = "N2 gene: WBGene00020062"), vjust = -1, color = '#3B2F2F', size = 5, fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_text(size = 20, color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(fill = NA)) +
    # plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "AB1 genome coordinates (Mb)", x = "N2 genome coordinates (Mb)")
syn6

ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/syntenicGene.png", dpi = 600, syn6, width = 14, height = 12)



syn2 <- joined %>% dplyr::filter(strain == "CB4856") %>%
  dplyr::arrange(start_gene) %>%
    dplyr::group_by(contig, n2_gene) %>%
    dplyr::filter(n() == 1) %>%
    dplyr::ungroup() 

  
syn7 <- ggplot(syn2 %>% dplyr::filter(contig == "ptg000026l" & (n2_gene == "WBGene00019185" | n2_gene == "WBGene00021043" | n2_gene == "WBGene00016222"))) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "#DB6333", alpha = 0.3) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.7) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6), color = "blue", linewidth = 1) +
  # geom_text(data = syn, 
            # aes(x =  (start_aln + 300) / 1e6 , y = (i.end  - (i.end - i.start + 225) / 2) / 1e6, label = paste0("WS gene: ", attributes)), vjust = -1, color = '#3B2F2F', size = 4.5, fontface = 'bold') +
  # geom_text(data = syn, 
            # aes(x =  (end_gene - (end_gene - start_gene) / 2) / 1e6 , y = WSS / 1e6, label = paste0("N2 gene: ", n2_gene)), vjust = -1, color = '#3B2F2F', size = 4.5, fontface = 'bold') +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(fill = NA)) +
  # plot.title = element_text(size = 14, color = 'black')) + 
  labs(y = "AB1 genome coordinates (Mb)", x = "N2 genome coordinates (Mb)")
syn7
# test <- joined %>% dplyr::filter((end_gene - start_gene) > 4e3) %>% dplyr::filter(strain == "CB4856")


ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/syntenicGene_three.png", dpi = 600, syn6, width = 14, height = 12)




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
# 






