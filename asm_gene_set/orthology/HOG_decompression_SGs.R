library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

syntelog_matrix <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/syntelog_matrix.tsv")


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
all_relations_rowid <- all_relations %>% dplyr::rename(rowid = HOG)
ortho_genes_dd_rowid <- all_ortho_genes_dd %>% dplyr::rename(rowid = HOG)

# Extract rowids for HOGs
complex_all_rowids <- all_relations_rowid %>% 
  dplyr::filter(if_any(2:ncol(.), ~ . > 1)) %>% 
  dplyr::pull(rowid) 

complex_rowids <- setdiff(complex_all_rowids, private_rowids) #remove complex HOGs that are PRIVATE!

complex_HOGS <- ortho_genes_dd_rowid %>% 
  dplyr::filter(rowid %in% complex_rowids) %>%
  dplyr::rename(HOG = rowid) 

# Splitting complex, PRIVATE HOGs into individual rows of genes
simple_HOGs <- ortho_genes_dd_rowid %>%
  dplyr::filter(!rowid %in% complex_rowids) %>%
  dplyr::rename(HOG = rowid) %>%
  tidyr::separate_rows(dplyr::everything(), sep = ',')


# ======================================================================================================================================================================================== #
# Decompress complex_HOGs based on syntelog matrix, and then dplyr::bind_rows() with simple_HOGs
# ======================================================================================================================================================================================== #


# a little test to get things working....
complex_subset <- complex_HOGS %>%
  dplyr::filter(!is.na(N2)) %>%
  dplyr::arrange(N2) %>%
  dplyr::slice_head(n=10) # grabbing some rows that have a single N2 gene in them 

N2_gene_subset <- complex_subset %>%
  dplyr::select(N2) %>%
  dplyr::pull()

syntelog_subset <- syntelog_matrix %>%
  dplyr::filter(n2_gene %in% N2_gene_subset) %>%
  dplyr::rename(N2 = n2_gene)

for N2_gene in syntelog_subset:
  when N2_gene == complex_subset$N2:
    for strain in N2_gene:
      when strain == 









