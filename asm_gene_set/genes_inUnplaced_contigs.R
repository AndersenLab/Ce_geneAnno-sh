library(dplyr)
library(ggplot2)
library(readr)

contig_breaks <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/henrique_lab/13_founders_contigsBreaks.tsv", col_names = c("contig","strain")) %>%
  tidyr::separate(contig, into = c("contig","break_start", "break_end"), sep = '_') %>%
  dplyr::mutate(break_start = as.numeric(break_start), break_end = as.numeric(break_end)) %>%
  dplyr::mutate(strain = ifelse(strain == "PB306", "ECA259", strain)) 

broken_ctgs <- contig_breaks %>% dplyr::distinct(contig,strain) %>% dplyr::mutate(broken = T)

contigs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/henrique_lab/13_founders_unplaced_contigs.tsv", col_names = c("contig","strain")) %>%
  dplyr::mutate(unplaced = T) %>% 
  dplyr::mutate(strain = ifelse(strain == "PB306", "ECA259", strain)) %>%
  dplyr::left_join(broken_ctgs, by = c("strain","contig")) %>% 
  dplyr::filter(is.na(broken)) %>%
  dplyr::select(-broken)

strains <- contigs %>% dplyr::distinct(strain) %>% dplyr::pull()
genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/142_140WSs_andCGC1_longestIsoGenes.tsv", col_names = c("contig","tool", "type","start","end","score","strand","phase","gene","strain"))

genes_in_unplaced <- genes %>%
  dplyr::left_join(contigs, by = c("contig","strain")) %>%
  dplyr::filter(!is.na(unplaced)) %>%
  dplyr::distinct(contig,start,end, gene,strain) 

genes_in_broken <- contig_breaks %>%
  dplyr::left_join(genes, by = c("contig","strain")) %>%
  dplyr::filter(start <= break_end & end >= break_start) %>%
  dplyr::distinct(contig,start,end,gene,strain)

final_genes_in_unplaced <- dplyr::bind_rows(genes_in_unplaced, genes_in_broken) %>% dplyr::arrange(strain)

contigs_per_strain <- final_genes_in_unplaced %>%
  dplyr::distinct(strain, contig) %>%
  dplyr::count(strain, name = "number_ctgs")

stats <- final_genes_in_unplaced %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(strain_gene_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(contigs_per_strain, by = "strain") %>%
  dplyr::distinct(strain, strain_gene_count, number_ctgs)


stats_long <- stats %>%
  # dplyr::mutate(strain_gene_count = strain_gene_count / 1000) %>%
  tidyr::pivot_longer(
    cols = c(strain_gene_count, number_ctgs),
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(metric = recode(metric, strain_gene_count = "Genes", number_ctgs = "Contigs"))

ggplot(stats_long, aes(x = strain, y = value, fill = metric)) +
  scale_fill_manual(values = c("Genes" = "forestgreen", "Contigs" = "blue")) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = value, group = metric), 
            position = position_dodge(width = 0.85), 
            vjust = -1, hjust = 0.5, size = 4, color = "black") +
  scale_y_continuous(name = "Number of predicted PC genes in unplaced contigs",
                     sec.axis = sec_axis(~ . / 1, name = "Number of unplaced contigs with predicted PC genes"), expand = c(0,1)) +
  labs(x = NULL, fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = 'black', face = 'bold', angle = 60, hjust = 1, size = 14),
    axis.text.y = element_text(color = 'black', size = 12),
    legend.text = element_text(color  = 'black', size = 18),
    axis.title.y = element_text(color = 'black', size = 16, face = 'bold'),
    legend.position = "top")

################## WHAT IS THE TOTAL AMOUNT OF SEQUENCE CONTAINED IN THESE UNPLACED CONTIGS??? ###########################
######### The number of genes is inflated because some of these contigs only PART of the contig is unplaced..... I need to filter for coordinates!!!!!!

class_genes <- final_genes_in_unplaced %>%
  dplyr::mutate(gene = sub("ID=","",gene)) %>%
  dplyr::mutate(in_unplaced = T) %>%
  dplyr::select(-contig)

strains <- contigs %>% dplyr::distinct(strain) %>% dplyr::pull()

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

all_class <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::left_join(class, by = "Orthogroup") %>%
  dplyr::select(Orthogroup, any_of(strains), class)

long_class <- all_class %>%
  tidyr::pivot_longer(
    cols = -c(Orthogroup, class),
    names_to = "strain",
    values_to = "gene",
    values_drop_na = TRUE) %>%
  tidyr::separate_rows(gene, sep = ",\\s*") %>%
  dplyr::select(strain, gene, class) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene))

prop <- class_genes %>% 
  dplyr::left_join(long_class, by = c("strain","gene")) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_genes = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(strain,class) %>%
  dplyr::mutate(num_in_class = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prop_unplaced_genes_class = num_in_class / num_genes * 100) %>%
  tidyr::complete(strain, class = c("core", "accessory", "private"),fill = list(num_in_class = 0, prop_unplaced_genes_class = 0)) %>%
  dplyr::distinct(strain,num_in_class,prop_unplaced_genes_class,class)

prop <- prop %>% mutate(class = factor(class, levels = c("core", "accessory", "private")))

ggplot(prop, aes(x = strain, y = prop_unplaced_genes_class, fill = class)) +
  geom_col(position = position_dodge(width = 0.85),width = 0.8) +
  geom_text(aes(label = num_in_class), position = position_dodge(width = 0.85), vjust = -1, hjust = 0.6, size = 4, color = "black") +
  scale_fill_manual(values = c(core = "green4",accessory = "#DB6333",private = "magenta3")) +
  scale_y_continuous(name = "Proportion of genes on unplaced contigs (%)",limits = c(0, 100),expand = c(0, 5)) +
  labs(x = NULL, fill = "Gene class") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 60, hjust = 1, face = "bold", size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14),
    legend.position = "top")


# Any single-copy orthologs in the core set lost???

sgo <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups_SingleCopyOrthologues.txt", col_names = "Orthogroup") %>%
  dplyr::pull(Orthogroup)

class2 <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(Orthogroup,class) %>%
  dplyr::filter(Orthogroup %in% sgo)

all_class2 <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::left_join(class2, by = "Orthogroup") %>%
  dplyr::select(Orthogroup, any_of(strains), class) %>%
  dplyr::filter(!is.na(class))

long_class2 <- all_class2 %>%
  tidyr::pivot_longer(
    cols = -c(Orthogroup, class),
    names_to = "strain",
    values_to = "gene",
    values_drop_na = TRUE) %>%
  tidyr::separate_rows(gene, sep = ",\\s*") %>%
  dplyr::select(strain, gene, class) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene))

prop2 <- class_genes %>% 
  dplyr::left_join(long_class2, by = c("strain","gene")) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_genes = n()) %>%
  dplyr::ungroup() 





