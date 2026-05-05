library(ggplot2)
library(dplyr)
library(readr)

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




