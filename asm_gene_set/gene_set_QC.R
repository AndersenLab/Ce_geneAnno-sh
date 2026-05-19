library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggrepel)


#######################################################################################################################################
### Looking at the coding sequence length of each gene set and proportion of single-exon genes in each strain
#######################################################################################################################################
# Loading in data
genes_sets <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv") %>%
  dplyr::select(-Orthogroup)

gff <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/140_CDS_features.tsv", col_names = c("type","start","end","strain","gene")) %>% dplyr::filter(strain != "ECA396") %>%
  dplyr::select(-type) 


# Calculating coding sequence length for all genes in each gene set
all_genes_gff <- gff %>%
  dplyr::mutate(span = end - start) %>%
  dplyr::group_by(gene,strain) %>%
  dplyr::mutate(gene_length = sum(span)) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_length, gene, strain) %>%
  dplyr::distinct()

all_genes_sets <- all_genes_gff %>%
  dplyr::left_join(genes_sets, by = c("strain","gene")) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(avg_gene_length = mean(gene_length)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(class = paste0(class, " (",round(avg_gene_length, 0),")"))

ggplot(all_genes_sets) + 
  geom_histogram(aes(x = gene_length, fill = class), alpha = 0.3, bins = 500) + 
  scale_fill_manual(values = c("core (1314)" = "green4", "accessory (935)" = "#DB6333", "private (1073)" = "magenta3")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 20, color = 'black'),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.text = element_text(size = 18, color = 'black'),
    legend.title = element_text(size = 18, color = 'black')) +
  labs(x = "Coding length", y = "Gene count", fill = "Gene set") +
  scale_y_log10(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))

# Looking at the smallest genes....
small <- all_genes_sets %>%
  dplyr::select(-avg_gene_length) %>%
  dplyr::mutate(class = ifelse(grepl("core",class),"core", 
                               ifelse(grepl('accessory', class), "accessory", "private"))) %>%
  dplyr::filter(gene_length <= 100) %>%
  dplyr::group_by(class) %>%
  dplyr::mutate(small_count = n()) %>%
  dplyr::distinct(class, small_count)

smallest <- all_genes_sets %>%
  dplyr::select(-avg_gene_length) %>%
  dplyr::mutate(class = ifelse(grepl("core",class),"core", 
                               ifelse(grepl('accessory', class), "accessory", "private"))) %>%
  dplyr::filter(gene_length <= 100) %>%
  dplyr::group_by(class) %>%
  dplyr::arrange(gene_length) %>%
  dplyr::slice_head(n = 20)


# Looking at the number of exons and number of genes for each gene set
all_genes_gff_exon_num <- gff %>%
  dplyr::left_join(genes_sets, by = c("strain","gene")) %>%
  # dplyr::mutate(span = end - start) %>%
  dplyr::group_by(gene,strain) %>%
  # dplyr::mutate(gene_length = sum(span)) %>%
  dplyr::mutate(num_exons = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,gene,class,num_exons) %>%
  dplyr::select(strain, class, num_exons) %>%
  dplyr::group_by(class, num_exons) %>%
  dplyr::mutate(num_genes_exon_num = n()) %>%
  dplyr::ungroup() 

final_exons <- all_genes_gff_exon_num %>%
  dplyr::distinct(class, num_exons, num_genes_exon_num) %>%
  dplyr::arrange(num_exons)

ggplot(final_exons) +
  geom_col(aes(x = num_exons, y = num_genes_exon_num, fill = class)) +
  scale_fill_manual(values = c("core" = "green4", "accessory" = "#DB6333", "private" = "magenta3")) +
  facet_wrap(~class, scales = "free") +
  theme(
    strip.text = element_text(size = 20, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 16, color = "black"),
    legend.position = "none",
    axis.title = element_text(size = 20, color = 'black')) +
  labs(x = "Number of exons", y = "Gene count") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))
  

# CDF
cdf_exons <- all_genes_gff_exon_num %>%
  distinct(class, num_exons, num_genes_exon_num) %>%
  arrange(class, num_exons) %>%
  group_by(class) %>%
  mutate(
    cum_genes = cumsum(num_genes_exon_num),
    total_genes = sum(num_genes_exon_num),
    cum_prop = cum_genes / total_genes
  ) %>%
  ungroup()

zero_rows <- final_exons %>%
  distinct(class) %>%
  mutate(
    num_exons = 0,
    num_genes_exon_num = 0
  )

final_exons2 <- bind_rows(final_exons, zero_rows) %>%
  arrange(class, num_exons)

cdf_exons <- final_exons2 %>%
  group_by(class) %>%
  arrange(num_exons, .by_group = TRUE) %>%
  mutate(
    cum_genes = cumsum(num_genes_exon_num),
    total_genes = sum(num_genes_exon_num),
    cum_prop = cum_genes / total_genes
  ) %>%
  ungroup()

ggplot(cdf_exons) +
  geom_step(aes(x = num_exons, y = cum_prop * 100, color = class), linewidth = 1.2) +
  scale_color_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  )) +
  # facet_wrap(~class) +
  theme(
    strip.text = element_text(size = 20, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 16, color = "black"),
    legend.position = "none",
    axis.title = element_text(size = 20, color = 'black')
  ) +
  labs(x = "Number of exons", y = "Cumulative proportion of genes (%)") +
  scale_x_continuous(breaks = seq(0,30,2), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  coord_cartesian(xlim = c(0,30))


# Number of genes and number of single-exon genes:
num_genes <- gff %>%
  dplyr::distinct(strain,gene) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(gene_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, gene_count)

gene_count_single_exon_count <- gff %>%
  dplyr::left_join(num_genes, by = "strain") %>%
  dplyr::group_by(gene, strain) %>%
  dplyr::mutate(num_exons = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(num_exons == 1) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_single_exon_genes = n()) %>%
  dplyr::distinct(strain, gene_count, gene, num_exons, num_single_exon_genes) %>%
  dplyr::select(strain, gene_count, num_exons, num_single_exon_genes) %>%
  dplyr::mutate(single_exon_removed = gene_count - num_single_exon_genes)

gene_exon_count_plt <- gene_count_single_exon_count %>%
  dplyr::distinct(strain, gene_count, num_single_exon_genes, single_exon_removed) %>%
  tidyr::pivot_longer(
    cols = c(gene_count, num_single_exon_genes, single_exon_removed),
    names_to = "count_type",
    values_to = "count") %>%
  dplyr::mutate(count_type = recode(count_type,"gene_count" = "Total genes", "num_single_exon_genes" = "Single exon genes", "single_exon_removed" = "Single exon removed")) 

# Add N2 !!!!!!!!!!!!!!!!!!!!!!!!
n2_CDSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.CDS.nuclear.cleaned.gff3", 
                           col_names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "gene", "strain")) %>%
  dplyr::select(type,gene,strain) %>% 
  dplyr::mutate(gene_count = length(unique(gene))) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(exon_num = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(exon_num == 1) %>%
  dplyr::mutate(num_single_exon_genes = length(unique(gene)),
                single_exon_removed = gene_count - num_single_exon_genes) %>%
  dplyr::select(strain, gene_count, num_single_exon_genes, single_exon_removed) %>%
  tidyr::pivot_longer(cols = c(gene_count, num_single_exon_genes, single_exon_removed),
                      names_to = "count_type",
                      values_to = "count") %>%
  dplyr::mutate(count_type = recode(count_type,"gene_count" = "Total genes", "num_single_exon_genes" = "Single exon genes", "single_exon_removed" = "Single exon removed")) %>%
  dplyr::arrange(desc(count))

gene_exon_count_plt_final <- gene_exon_count_plt %>%
  dplyr::bind_rows(n2_CDSs) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(count_type = factor(count_type, levels = c("Total genes", "Single exon removed","Single exon genes"))) 

strain_order <- gene_exon_count_plt_final %>%
  dplyr::filter(count_type == "Total genes") %>%
  dplyr::distinct(strain, count) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::pull(strain)

gene_exon_count_plt_final <- gene_exon_count_plt_final %>%
  dplyr::mutate(strain = factor(strain, levels = strain_order))


ggplot(data = gene_exon_count_plt_final, aes(x = strain, y = count)) +
  geom_col(aes(fill = count_type), position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(aes(color = count_type), position = position_dodge(width = 0.8), size = 2, show.legend = FALSE) +
  geom_line(aes(color = count_type, group = count_type), position = position_dodge(width = 0.8), linewidth = 2, show.legend = FALSE) +
  labs(y = "Gene count", fill = "") +
  theme(
    strip.text = element_text(size = 20, color = 'black'),
    panel.background = element_blank(),
    legend.title = element_text(size = 18, color = 'black'),
    legend.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 20, color = 'black'),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 13),
    legend.position = "top") +
  scale_fill_manual(values = c("Total genes" = "black", "Single exon genes" = "blue", "Single exon removed" = "red")) +
  scale_color_manual(values = c("Total genes" = "black", "Single exon genes" = "blue", "Single exon removed" = "red")) +
  scale_y_continuous(breaks = seq(0,22000, 2000), expand = c(0,100))
















#######################################################################################################################################
# Looking at how many non-N2 genes do WSs have on average????
#######################################################################################################################################
# ======================================================================================================================================================================================== #
# Pulling all genes, coordinates, and alignments for all WSs and N2 #
# ======================================================================================================================================================================================== #
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/142_140WSs_andCGC1_longestIsoGenes.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
genes_strain <- rbind(genes_strain,N2_gff)
all_genes_strain <- genes_strain %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene")

N2_tranGene <- N2_gff %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(attributes = gsub("ID=transcript:","",attributes), attributes = gsub("Parent=gene:","",attributes)) %>%
  tidyr::separate_wider_delim(attributes, delim = ";",names = c("tran", "gene", "rest"), too_many = "merge") %>%
  dplyr::select(tran,gene, -rest)

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-IDY,-LENR,-LENQ) %>% dplyr::filter(strain != "ECA396")


# ======================================================================================================================================================================================== #
# HOG matrix manipulation and plotting #
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
  dplyr::count(freq, sum, class) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) 



gs_allOrtho <- ggplot(data = classification, aes(x = sum, y = n, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("Orthogroups") + # longest isoform
  xlab("Genomes") +
  # scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  # scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(classification$sum), by = 25)) +
  theme(
    axis.title = element_text(size = 24, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    plot.margin = margin(l = 20, r = 20, t = 20),
    # plot.title = element_text(size=26, face = 'bold', hjust=0.5),
    legend.text = element_text(size=22, color = 'black'),
    axis.text = element_text(size=18, color = 'black')
  )
gs_allOrtho

OG_class_count <- classification %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_OG = sum(n)) %>%
  dplyr::ungroup()


#######################################################################################################################################
### Looking at non-N2 ###
#######################################################################################################################################
all_relations_clean <- all_relations %>% dplyr::rename_with(~ gsub("_count", "", .x))

names <- colnames(all_relations_clean %>% dplyr::select(-N2, -Orthogroup))

nonRefGenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(nonRefGenes) <- c("strain","nonN2_genes",'nonN2_Orthogroups')

OG_list <- list()

for (i in 1:length(names)) {
  soi <- names[i]
  print(paste0("On strain: ", soi, ". ", i, "/141."))
  
  nonRefGenes[i,1] = soi
  
  non_ref <- all_relations_clean %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(is.na(N2) & !is.na(.data[[soi]])) %>% dplyr::select(.data[[soi]]) %>% sum(na.rm = TRUE) 
  
  nonRefGenes[i,2] = non_ref
  
  non_ref_og_count <- all_relations_clean %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(is.na(N2) & !is.na(.data[[soi]])) %>% nrow()
  non_ref_og <- all_relations_clean %>% dplyr::select(Orthogroup, N2, .data[[soi]]) %>% dplyr::filter(is.na(N2) & !is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
  OG_list[[i]] <- non_ref_og
  
  nonRefGenes[i,3] = non_ref_og_count
  
}

nonRefGenes_long_n2 <- nonRefGenes %>%
  pivot_longer(
    cols = c(nonN2_genes, nonN2_Orthogroups),
    names_to = "metric",
    values_to = "count"
  )

label_df_top_n2 <- nonRefGenes_long_n2 %>%
  dplyr::filter(metric == "nonN2_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)

label_df_bottom_n2 <- nonRefGenes_long_n2 %>%
  dplyr::filter(metric == "nonN2_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)

labels_df_n2 <- label_df_top_n2 %>% dplyr::bind_rows(label_df_bottom_n2)

n2 <- ggplot(nonRefGenes_long_n2, aes(x = metric, y = count)) +
  geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df_n2, aes(label = strain), size = 4, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("nonN2_genes" = "gray70", "nonN2_Orthogroups" = "gray30")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black')) +
  coord_cartesian(ylim = c(800,5500))


# Looking at the proportion in each gene set
OG_vector_N2 <- unique(unlist(OG_list))

OG_classes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv") %>% dplyr::select(Orthogroup, class) %>% dplyr::distinct() %>%
  dplyr::filter(Orthogroup %in% OG_vector_N2) 

OG_propN2 <- OG_classes %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "private"))) %>%
  dplyr::mutate(source = "N2")

ggplot(OG_prop, aes(x = prop, y = "", fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black'))


#######################################################################################################################################
### Looking at non-CGC1 ###
#######################################################################################################################################
all_relations_clean <- all_relations %>% dplyr::rename_with(~ gsub("_count", "", .x))

names <- colnames(all_relations_clean %>% dplyr::select(-CGC1,-Orthogroup))

nonRefGenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(nonRefGenes) <- c("strain","nonCGC1_genes",'nonCGC1_Orthogroups')

OG_list <- list()

for (i in 1:length(names)) {
  soi <- names[i]
  print(paste0("On strain: ", soi, ". ", i, "/141."))
  
  nonRefGenes[i,1] = soi
  
  non_ref <- all_relations_clean %>% dplyr::select(CGC1, .data[[soi]]) %>% dplyr::filter(is.na(CGC1) & !is.na(.data[[soi]])) %>% dplyr::select(.data[[soi]]) %>% sum(na.rm = TRUE) 
  
  nonRefGenes[i,2] = non_ref
  
  non_ref_og_count <- all_relations_clean %>% dplyr::select(CGC1, .data[[soi]]) %>% dplyr::filter(is.na(CGC1) & !is.na(.data[[soi]])) %>% nrow()
  non_ref_og <- all_relations_clean %>% dplyr::select(Orthogroup, CGC1, .data[[soi]]) %>% dplyr::filter(is.na(CGC1) & !is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
  OG_list[[i]] <- non_ref_og
  
  nonRefGenes[i,3] = non_ref_og_count
  
}

nonRefGenes_long <- nonRefGenes %>%
  pivot_longer(
    cols = c(nonCGC1_genes, nonCGC1_Orthogroups),
    names_to = "metric",
    values_to = "count"
  )

label_df_top <- nonRefGenes_long %>%
  dplyr::filter(metric == "nonCGC1_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)

label_df_bottom <- nonRefGenes_long %>%
  dplyr::filter(metric == "nonCGC1_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)

labels_df <- label_df_top %>% dplyr::bind_rows(label_df_bottom)

cg <- ggplot(nonRefGenes_long, aes(x = metric, y = count)) +
  geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df, aes(label = strain), size = 4, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("nonCGC1_genes" = "gray70", "nonCGC1_Orthogroups" = "gray30")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    # axis.title.y = element_text(size = 16, color = 'black'),
    # axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black')) +
  coord_cartesian(ylim = c(800,5500))


n2_cg <- cowplot::plot_grid(
  n2, cg,
  nrow = 1)
n2_cg

# Looking at the proportion in each gene set?
OG_vector_CGC1 <- unique(unlist(OG_list))

OG_classes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv") %>% dplyr::select(Orthogroup, class) %>% dplyr::distinct() %>%
  dplyr::filter(Orthogroup %in% OG_vector_CGC1) 

OG_prop <- OG_classes %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "private"))) %>%
  dplyr::mutate(source = "CGC1")

ggplot(OG_prop, aes(x = prop, y = "", fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black'))


both <- OG_prop %>% dplyr::bind_rows(OG_propN2)

both_prop <- ggplot(both, aes(x = prop, y = source, fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black'))


all_plt <- cowplot::plot_grid(
  n2_cg, both_prop,
  nrow = 2,
  rel_heights = c(4,1))
all_plt











#######################################################################################################################################
### Looking at N2-specific ###
#######################################################################################################################################
all_relations_clean <- all_relations %>% dplyr::rename_with(~ gsub("_count", "", .x))

names_all <- colnames(all_relations_clean %>% dplyr::select(-N2, -Orthogroup))

nonRefGenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(nonRefGenes) <- c("strain","N2_specific_genes",'N2_specific_Orthogroups')

OG_list <- list()

for (i in 1:length(names_all)) {
  soi <- names_all[i]
  print(paste0("On strain: ", soi, ". ", i, "/141."))
  
  nonRefGenes[i,1] = soi
  
  non_ref <- all_relations_clean %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% dplyr::select(N2) %>% sum(na.rm = TRUE) 
  
  nonRefGenes[i,2] = non_ref
  
  non_ref_og_count <- all_relations_clean %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% nrow()
  non_ref_og <- all_relations_clean %>% dplyr::select(Orthogroup, N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
  OG_list[[i]] <- non_ref_og
  
  nonRefGenes[i,3] = non_ref_og_count
  
}

nonRefGenes_long <- nonRefGenes %>%
  pivot_longer(
    cols = c(N2_specific_genes, N2_specific_Orthogroups),
    names_to = "metric",
    values_to = "count"
  )

label_df_top <- nonRefGenes_long %>%
  dplyr::filter(metric == "N2_specific_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)

label_df_bottom <- nonRefGenes_long %>%
  dplyr::filter(metric == "N2_specific_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)

labels_df <- label_df_top %>% dplyr::bind_rows(label_df_bottom)

n2_spec <- ggplot(nonRefGenes_long, aes(x = metric, y = count)) +
  geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
  geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df, aes(label = strain), size = 4, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("N2_specific_genes" = "orange", "N2_specific_Orthogroups" = "#DB6333")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16, color = 'black')) #+
  # coord_cartesian(ylim = c(800,5500))
n2_spec

# Looking at proprotion in each gene set
OG_vector_N2_spec <- unique(unlist(OG_list))

OG_classes_spec <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv") %>% dplyr::select(Orthogroup, class) %>% dplyr::distinct() %>%
  dplyr::filter(Orthogroup %in% OG_vector_N2_spec) 

OG_propN2_spec <- OG_classes_spec %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "private"))) %>%
  dplyr::mutate(source = "N2")

n2_spec_prop <- ggplot(OG_propN2_spec, aes(x = prop, y = "", fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black')) #+
  # scale_x_continuous(expand = c(0,0)) +
  # scale_y_discrete(expand = c(0,0))

n2_specific_genes <- cowplot::plot_grid(
  n2_spec, n2_spec_prop,
  nrow = 2,
  rel_heights = c(6,1))
n2_specific_genes


