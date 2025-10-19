library(plyr) # ALWAYS LOAD BEFORE DPLYR
library(readr)
library(org.Ce.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(ape)
library(tidyr)
library(data.table)
library(tibble)
library(clusterProfiler) ## BiocManager::install("clusterProfiler") # need this and the next package???
library(enrichplot)
library(cowplot)
library(mgcv)


# ======================================================================================================================================================================================== #
# Creating a key-value for transcripts and gene names among all strain of C. elegans (115 WS and N2) for converting transcripts to genes in the N0.tsv output from Orthofinder #
# ======================================================================================================================================================================================== #

# gff <- ape::read.gff("/vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3")
# 
# # Load in concatenated GFF of N2 and all WSs
# all_gff <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/0507_116_geneAndmRNA_gff.tsv", col_names = c("type", "start", "end", "strand", "attributes"))
# 
# all_gff_formatted <- all_gff %>%
#   dplyr::mutate(attributes = gsub("ID=","", attributes)) %>%
#   dplyr::mutate(attributes_new = str_extract(attributes, "transcript:[^;]+"), parent = str_extract(attributes, "Parent=gene:[^;]+")
#   ) %>%
#   dplyr::mutate(attributes_new = ifelse(is.na(attributes_new),attributes,attributes_new)) %>%
#   dplyr::mutate(parent = ifelse(is.na(parent),attributes,parent)) %>%
#   dplyr::select(-attributes) %>%
#   dplyr::mutate(parent = gsub("Parent=gene:","",parent)) %>%
#   dplyr::mutate(attributes_new = gsub("transcript:","",attributes_new)) %>%
#   dplyr::mutate(attributes_new = sub(";P.*", "", attributes_new), parent = sub(".*Parent=([^;]+);?.*", "\\1", parent)) %>%
#   dplyr::mutate(parent = sub(".*gene:(WBGene[0-9]+).*", "\\1", parent), attributes_new = sub(".*gene:(WBGene[0-9]+).*", "\\1", attributes_new)) %>%
#   dplyr::mutate(parent = gsub(";","",parent), attributes_new = gsub(";","",attributes_new))
# 
# transcript_key <- all_gff_formatted %>%
#   dplyr::filter(type == "mRNA") %>%
#   dplyr::select(attributes_new,parent) %>%
#   dplyr::rename(transcript = attributes_new, gene = parent)

# write.table(transcript_key,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/transcripts_gene_115WI_N2.tsv", quote = F, row.names = F, col.names = F, sep = '\t')





# ======================================================================================================================================================================================== #
# Pulling all genes, coordinates, and alignments for all WSs and N2 #
# ======================================================================================================================================================================================== #
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
all_genes_strain <- genes_strain %>%
  dplyr::filter(strain != "N2" | grepl("protein_coding", attributes)) %>% #filter out non-protein coding N2 genes
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) 

# test1 <- all_genes_strain %>%
#   dplyr::filter(strain == "N2") # ~19,000

# write.table(all_genes_strain,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/116strain_genes.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/nucmer_runs/115_WI_transformed.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-IDY,-LENR,-LENQ)
# write.table(nucmer,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/115WS_nucmer_clean.tsv", quote = F, row.names = F, col.names = T, sep = '\t')









##### Verifying that transcript to gene conversion works correctly ############################################################################################################
# orthogroups_tran <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/prot_115/OrthoFinder/Results_Apr22/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
# strCol <- colnames(orthogroups_tran)
# strCol_c1 <- gsub(".braker.protein","",strCol)
# strCol_c2 <- gsub("_WS283.protein","",strCol_c1)
# colnames(orthogroups_tran) <- strCol_c2
# 
# # print(nrow(orthogroups_tran)) # 64377
# 
# #### Counting number of transcripts per orthogroup for each strain #### 
# ortho_count_tran <- orthogroups_tran
# 
# strCol_c2_u <- strCol_c2[!strCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]
# 
# for (i in 1:length(strCol_c2_u)) {
#   print(paste0(i,"out of", length(strCol_c2_u)))
#   temp_colname = paste0(strCol_c2_u[i], "_count")
# 
#   ortho_count_tran <- ortho_count_tran %>%
#     dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strCol_c2_u[i]),", ") + 1)
# }
# 
# 
# all_relations_tran <- ortho_count_tran %>%
#   dplyr::select(HOG, dplyr::contains("_count"))
# 
# 
# #### Transcripts converted to genes #### see script: /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/tran_gene.sh
# ortho_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0422_genes.tsv")
# 
# stCol <- colnames(ortho_genes)
# stCol_c1 <- gsub(".braker.protein","",stCol)
# stCol_c2 <- gsub("_WS283.protein","",stCol_c1)
# colnames(ortho_genes) <- stCol_c2
# 
# ortho_count_dup <- ortho_genes
# 
# stCol_c2_u <- stCol_c2[!stCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]
# 
# for (i in 1:length(stCol_c2_u)) {
#   print(paste0(i,"out of", length(stCol_c2_u)))
#   temp_colname = paste0(stCol_c2_u[i], "_count")
#   
#   ortho_count_dup <- ortho_count_dup %>%
#     dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(stCol_c2_u[i]),", ") + 1)
# }
# 
# all_relations_duplicated_genes <- ortho_count_dup %>%
#   dplyr::select(HOG, dplyr::contains("_count"))
# 
# ## all_relations and all_relations_g should be identical
# identical(all_relations_tran, all_relations_duplicated_genes) # TRUE
#################################################################################################################################################################################








# 
# 
# N2_braker <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_N2_braker_genes.tsv")
# 
# 
# strainCol <- colnames(N2_braker)
# strainCol_c1 <- gsub(".braker.longest.protein","",strainCol)
# strainCol_c2 <- gsub(".longest.protein","",strainCol_c1)
# colnames(N2_braker) <- strainCol_c2
# 
# ortho_count <- N2_braker
# 
# strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]
# 
# for (i in 1:length(strainCol_c2_u)) {
#   print(paste0(i,"out of", length(strainCol_c2_u)))
#   temp_colname = paste0(strainCol_c2_u[i], "_count")
#   
#   ortho_count <- ortho_count %>%
#     dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
# }
# 
# all_relations <- ortho_count %>%
#   dplyr::select(HOG, dplyr::contains("_count"))
# 
# 
# #### Plotting classification based on all HOGs ####
# private_freq = (0.5)
# 
# classification <- all_relations %>%
#   dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
#   dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
#   dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
#   dplyr::mutate(
#     class = case_when(
#       freq == 1 ~ "core",
#       freq == private_freq ~ "private",
#       TRUE ~ "undefined"
#     )
#   ) %>%
#   dplyr::count(freq, class) %>%
#   dplyr::mutate(percent = (n / sum(n)) * 100) 
# 
# gs_allOrtho <- ggplot(data = classification, aes(x = freq * 100, y = percent, fill = class)) + 
#   geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
#   scale_fill_manual(values = c(
#     "core" = "green4",
#     "private" = "magenta3"
#   ), 
#   limits = c("core", "private"),  # Manually ordering legend items
#   guide = guide_legend(title = NULL) 
#   ) +
#   ylab("Percent of HOGs - Longest Isoform") +
#   xlab("Frequency") +
#   # ggtitle("All orthogroups") +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
#   scale_x_continuous(labels = scales::percent_format(scale = 1)) +
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 16),
#     legend.position = c(0.85, 0.8),
#     plot.title = element_text(size=18, face = 'bold', hjust=0.5),
#     legend.text = element_text(size=13, color = 'black'),
#     axis.text = element_text(size=12, color = 'black')
#   )
# gs_allOrtho
# 
# HOG_class_count <- classification %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarise(n_HOG = sum(n)) %>%
#   dplyr::ungroup()
# 
# N2_braker <- all_relations %>%
#   dplyr::select(N2_count) %>%
#   dplyr::filter(is.na(N2_count)) # 25 - number of private HOGs belonging to the N2.WS283
# 
# N2_WS283 <- all_relations %>%
#   dplyr::select(N2.WS283_count) %>%
#   dplyr::filter(is.na(N2.WS283_count)) # 96 - number of private HOGs belonging to the BRAKER predictions
# 
# n2_braker_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/N2_BRAKER_error/N2.BRAKER.longest.tsv", col_names = c("contig", "type", "start", "end", "strand", "attributes")) 
# ws283_genes <- all_genes_strain %>%
#   dplyr::filter(strain == "N2")
# 
# n2_braker_WS283_genes <- n2_braker_genes %>%
#   dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
#   dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
#   dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
#   dplyr::mutate(strain = "N2_BRAKER") %>%
#   dplyr::bind_rows(ws283_genes) 
# 
# n2_gene_count <- n2_braker_WS283_genes %>%
#   dplyr::count(strain, name = "n_genes")
# 
# ggplot(n2_gene_count) + 
#   geom_bar(aes(x = strain, y = n_genes, fill = strain), stat = "identity") +
#   geom_text(aes(x = strain, y = n_genes, label = n_genes), vjust = -0.3, size = 5) +
#   ylab("Number of PC genes") +
#   theme(
#     legend.position = 'none',
#     panel.background = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 14, face = 'bold'),
#     axis.text.y = element_text(size = 13),
#     axis.text.x = element_text(face = 'bold', size=14, color = 'black'))








# ======================================================================================================================================================================================== #
# HOG matrix manipulation and plotting #
# ======================================================================================================================================================================================== #

# Converting transcripts to genes and removing all duplicate genes (not needed anymore, using longest isoform) in a dataframe cell:
# see script - /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/tran_gene.sh 
# ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/N0_0504_genes.tsv")
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Oct16/Orthogroups/Orthogroups.tsv")


# whatever <- ortho_genes_dd %>%
  # dplyr::select(JU1581.braker.longest.protein, N2.longest.protein)

strainCol <- colnames(ortho_genes_dd)
ugh <- gsub(".20251012.inbred.blobFiltered.softMasked.braker.longestIso.protein","", strainCol)
ugh2 <- gsub("..20251014.inbred.blobFiltered.softMasked.braker.longestIso.protein","",ugh)
strainCol_c2 <- gsub("c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein","N2", ugh2)
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


private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Oct16/Orthogroups/Orthogroups_UnassignedGenes.tsv") 

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


## Create pangenome size and core set curve after iterations ##

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/gene_set_allOrtho_115.png", gs_allOrtho, height = 5, width = 11, dpi = 600)


### PLOTTING BASED ON GENES ###
classification_genes <- all_relations %>%
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
  ) 

sum_genes <- all_relations %>%
  dplyr::mutate(all_genes = rowSums(dplyr::across(2:ncol(.)), na.rm = TRUE)) %>%
  dplyr::select(all_genes) %>%
  dplyr::bind_cols(classification_genes) %>%
  dplyr::group_by(freq, sum, class) %>%
  dplyr::summarise(total_genes = sum(all_genes, na.rm = TRUE)) %>%
  dplyr::ungroup()


genes_allHOGs <- ggplot(data = sum_genes, aes(x = sum, y = total_genes / 1000, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("Genes (1e3)") + # longest isoform
  xlab("Strains") +
  # scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  # scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    # plot.margin = margin(l = 20),
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=16, color = 'black'),
    axis.text = element_text(size=14, color = 'black')
  )
genes_allHOGs

all_genes_class_count <- sum_genes %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_genes = sum(total_genes)) %>%
  dplyr::ungroup()





### PLOTTING BAR PLOTS FOR THE PROPORTION OF GENES THAT ARE CLASSIFIED AS EACH GENE SET IN STRAIN
table <- all_relations %>%
  dplyr::left_join(classification_genes %>% dplyr::select(Orthogroup, class), by = "Orthogroup") %>%
  dplyr::select(-Orthogroup)

colnames(table) <- gsub("_count", "", colnames(table))

count_cols <- colnames(table)[1:(ncol(table) - 1)] 

results_list <- list()

for (col_name in count_cols) {
  temp <- table %>%
    dplyr::select(dplyr::all_of(col_name), class) %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(n_genes = sum(!!sym(col_name), na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = class, values_from = n_genes)
  
  # Add column to identify source column
  temp <- temp %>%
    dplyr::mutate(strain = col_name) %>%
    dplyr::select(strain, core, accessory, private)  
  
  results_list[[col_name]] <- temp
}

final_df <- bind_rows(results_list)

df_long <- final_df %>%
  tidyr::pivot_longer(cols = c(core, accessory, private), names_to = "class", values_to = "n_genes") %>%
  dplyr::mutate(class = factor(class, levels = c("private", "accessory", "core")), strain = factor(strain, levels = rev(unique(strain))))

ugh <- df_long %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(meann = mean(n_genes))
# 
# core <- ugh %>%dplyr::filter(class == "core")
# acc <- ugh %>% dplyr::filter(class == "accessory") %>% dplyr::mutate(meann = meann + core$meann)
# priv <- ugh %>% dplyr::filter(class == "private") %>% dplyr::mutate(meann = (meann + acc$meann))
# 
# contrib <- ggplot(df_long, aes(x = n_genes, y = strain, fill = class), alpha = 0.4) +
#   scale_fill_manual(values = c(
#     core = "green4",
#     accessory = "#DB6333",
#     private = "magenta3"
#   )) +
#   geom_bar(stat = "identity") +
#   geom_vline(xintercept = core$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
#   geom_vline(xintercept = acc$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
#   geom_vline(xintercept = priv$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
#   labs(x = "Number of Genes", fill = "Gene set") +
#   theme(
#     axis.text.y = element_text(size = 8.5, color = 'black'),
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(size = 12, color = 'black'),
#     axis.title.x = element_text(size = 14, color = 'black', face = 'bold'),
#     panel.background = element_blank(),
#     panel.border = element_rect(color = 'black', fill = NA)) +
#   scale_x_continuous(expand = c(0, 250))
# contrib

df_percent <- final_df %>%
  dplyr::mutate(total = core + accessory + private) %>%
  dplyr::mutate(core = 100 * core / total, accessory = 100 * accessory / total, private = 100 * private / total) %>%
  dplyr::select(-total) %>%
  tidyr::pivot_longer(cols = c(core, accessory, private), names_to = "class", values_to = "percent")


strain_order <- df_percent %>%
  dplyr::filter(class == "core") %>%
  dplyr::arrange(percent) %>%
  dplyr::pull(strain)

df_percent <- df_percent %>%
  dplyr::mutate(class = factor(class, levels = c("private", "accessory", "core")), strain = factor(strain, levels = strain_order))

core <- ugh %>% dplyr::filter(class == "core")
acc <- ugh %>% dplyr::filter(class == "accessory") %>% dplyr::mutate(meann = meann + core$meann)
priv <- ugh %>% dplyr::filter(class == "private") %>% dplyr::mutate(meann = 100 - meann)


contrib <- ggplot(df_percent, aes(x = percent, y = strain, fill = class)) +
  scale_fill_manual(values = c(
    core = "green4",
    accessory = "#DB6333",
    private = "magenta3"
  )) +
  geom_bar(stat = "identity") +
  # geom_vline(xintercept = core$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
  # geom_vline(xintercept = acc$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
  # geom_vline(xintercept = priv$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
  labs(x = "Percent of genes", fill = "Gene set") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    axis.text.y = element_text(size = 8.5, color = 'black'),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 22, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)
  )
contrib

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/bar.png", contrib, height = 13, width = 12, dpi = 600)






# ======================================================================================================================================================================================== #
# Plotting N2 genes #
# ======================================================================================================================================================================================== #
# Turning gene count table into a binary presence/absense of an orthogroup 
count <- all_relations %>%
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
  ) 


n2_gene <- all_genes_strain %>%
  dplyr::filter(strain == "N2") %>%
  dplyr::rename(ID = attributes)

n2_table <- ortho_genes_dd %>% 
  dplyr::select(-OG,-'Gene Tree Parent Clade') %>%
  dplyr::bind_rows((private_OGs %>% dplyr::rename(HOG = Orthogroup))) %>%
  dplyr::select(HOG,N2) 

ortho_count_wCoord <- count %>%
  dplyr::left_join(n2_table, by = "HOG") %>%
  dplyr::select(freq, class, N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ",\\s*") %>% # splitting rows so each gene is on a row and it retains is gene set classification
  dplyr::left_join(n2_gene, by = c("N2" = "ID")) %>%
  dplyr::select(freq,class,N2,contig,start,end) 

HDRs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/114HDRs.tsv", col_names = c("CHROM","START","END","strain"))

all_regions <- HDRs %>%
  dplyr::arrange(CHROM,START) %>%
  dplyr::group_split(CHROM)

strain_count <- HDRs %>% dplyr::distinct(strain, .keep_all = T)
print(nrow(strain_count))


getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>% 
        dplyr::arrange(CHROM,START) %>%
        dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      #print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
      if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
        print("NO MORE INTERSECTS")
        k=0
      } else {
        
        temp <- checkIntersect %>%
          dplyr::mutate(gid=data.table::rleid(check)) %>%
          dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))
        
        collapse <- temp %>%
          dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart=min(START)) %>%
          dplyr::mutate(newEnd=max(END)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(START=newStart,END=newEnd) %>%
          dplyr::select(-newEnd,-newStart)
        
        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

HDR_collapse_master <- getRegFreq(all_regions)

all_collapsed <- ldply(HDR_collapse_master, data.frame) %>%
  dplyr::select(-strain) 


# Adding resolution if genes are found in a HDR or not
ortho_count_wCoord_HDR <- ortho_count_wCoord %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(in_HDR = any(
    contig == all_collapsed$CHROM &
      start >= all_collapsed$START &
      end <= all_collapsed$END
  )) 


stats <- ortho_count_wCoord_HDR %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::ungroup()

ortho_count_wCoord_HDR_final <- ortho_count_wCoord_HDR %>%
  dplyr::left_join(stats, by = 'class') %>%
  dplyr::filter(contig != "MtDNA") %>%
  dplyr::mutate(class = dplyr::recode(class, "core" = "core (14,799)", "accessory" = "accessory (4,660)", "private" = "private (462)")) %>%
  dplyr::rename(Class = class, HDR = in_HDR)


N2_genes_plot <- ggplot(ortho_count_wCoord_HDR_final) +
  geom_point(aes(x = start / 1e6, y = freq * 100, color = Class, shape = HDR), size = 2) +
  scale_color_manual(
    name = "Gene set",
    values = c("core (14,799)" = "green4", "accessory (4,660)" = "#DB6333", "private (462)" = "magenta3"),
    limits = c("core (14,799)", "accessory (4,660)", "private (462)")) +
  scale_shape_manual(
    name = "In a HDR?",
    values = c("TRUE" = 4, "FALSE" = 1)) +
  geom_smooth(aes(x = start / 1e6, y = freq * 100), method = "loess", se = TRUE, color = "blue") +
  facet_wrap(~contig, ncol = 6, scales = "free_x") +
  labs(x = "Genome position (Mb)", y = "N2 gene frequency") +
  scale_y_continuous(expand = c(0.01, 0), labels = scales::percent_format(scale = 1)) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 16, color = "black")
  )
N2_genes_plot

core_hdr <- ortho_count_wCoord_HDR_final %>%
  dplyr::filter(Class == "core (14,799)") %>%
  dplyr::count(HDR) # 1958/12841:InHDRs/notInHDRs - 13%

acc_hdr <- ortho_count_wCoord_HDR_final %>%
  dplyr::filter(Class == "accessory (4,660)") %>%
  dplyr::count(HDR) # 2273/2380:InHDRs/notInHDRs - 49%

priv_hdr <- ortho_count_wCoord_HDR_final %>%
  dplyr::filter(Class == "private (462)") %>%
  dplyr::count(HDR) # 105/357:InHDRs/notInHDRs - 23%

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/N2_geneFreq.png", N2_genes_plot, height = 10, width = 15, dpi = 600)


ortho_count_wCoord_HDR_final$HDR <- factor(ortho_count_wCoord_HDR_final$HDR, levels = c(TRUE, FALSE))

ggplot(ortho_count_wCoord_HDR_final) +
  geom_jitter(aes(x = HDR, y = freq * 100), width = 0.3, alpha = 0.5, size = 1) +
  geom_boxplot(aes(x = HDR, y = freq * 100, fill = HDR), outlier.alpha = 0, alpha = 0.5) +
  scale_fill_manual(
    name = "Hyper-divergent region?",  
    labels = c("non-HDR", "HDR"),  
    values = c("#6495ED","#FF6347")  
  ) +
  labs(x="", y = "Frequency of all N2 genes in a HOG") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        panel.grid = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))
        # legend.key.height = unit(0.3, "cm"),  
        # legend.key.width = unit(0.3, "cm"))


 # SOME TYPE OF VISUALIZATION OF GENES FOUND IN N2 VERSUS NOT IN N2??








# ======================================================================================================================================================================================== #
# Plotting rarefaction #
# ======================================================================================================================================================================================== #
# private <- all_relations %>%
#   dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
#   dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
#   dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
#   dplyr::mutate(
#     class = case_when(
#       freq == 1 ~ "core",
#       freq >= 0.95 & freq < 1 ~ "soft-core",
#       freq > private_freq & freq < 0.95 ~ "accessory",
#       freq == private_freq ~ "private",
#       TRUE ~ "undefined"
#     )
#   ) %>%
#   dplyr::filter(class == "private")
# 
# strains <- private %>%
#   dplyr::select(-HOG,-sum,-freq,-class) %>%
#   colnames()
# 
# private_ordered <- private %>%
#   dplyr::select(all_of(strains)) %>%  # Keep only strain columns
#   dplyr::summarise(across(everything(), sum, na.rm = TRUE)) %>%
#   pivot_longer(cols = everything(), names_to = "strain", values_to = "count") %>%  # Use a string for `names_to`
#   dplyr::arrange(desc(count)) %>%
#   dplyr::pull(strain)
# 
# private_final <- private %>%
#   dplyr::select(HOG, all_of(private_ordered), sum, freq, class)
# 
# rarefaction <- data.frame(
#   num_strains = integer(),
#   num_private_orthogroups = integer()
# )
# 
# for (i in 2:length(strainCol_c2_u)) {
#    temp <- private_final %>% dplyr::select(2:i)
#    # print(temp)
# 
#    private_count <- sum(apply(temp, 1, function(row) sum(!is.na(row)) == 1))
#    rarefaction <- rbind(rarefaction, data.frame(num_strains = i-1, num_private_orthogroups = private_count))
# }
# 
# rfc <- ggplot(data = rarefaction) +
#   geom_point(aes(x=num_strains, y = num_private_orthogroups), size=2.5, color = 'magenta3', alpha=0.8) +
#   xlab("Genomes") +
#   ylab("Private HOGs") +
#   # coord_cartesian(ylim = c(200,1500)) +
#   theme(
#     panel.background = element_blank(),
#     axis.title = element_text(size=16, color = 'black', face = 'bold'),
#     axis.text = element_text(size=14, color = 'black'),
#     panel.border = element_rect(fill = NA),
#     plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  # 20pt on all sides
#   )
# 
# rfc
# 
# 
# # ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/rarefaction_115_labelled.png", height = 5, width = 9, rfc_labeled, dpi = 600)


### PLOTTING RAREFRACTION BASED ON PRIVATE GENE NUMBERS NOT PRIVATE HOGs ###
private <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq >= 0.95 & freq < 1 ~ "soft-core",
      freq > private_freq & freq < 0.95 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>%
  dplyr::filter(class == "private")

priv_OGs <- private %>%
  dplyr::pull(Orthogroup)

strains <- private %>%
  dplyr::select(-Orthogroup,-sum,-freq,-class) %>%
  colnames()

priv_genes <- all_relations %>%
  dplyr::filter(Orthogroup %in% priv_OGs) 

private_ordered_genes <- priv_genes %>%
  dplyr::select(all_of(strains)) %>%  # Keep only strain columns
  dplyr::summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "strain", values_to = "count") %>% 
  dplyr::arrange(desc(count)) #%>%
  dplyr::pull(strain)

private_final <- priv_genes %>%
  dplyr::select(Orthogroup, all_of(private_ordered_genes))

rarefaction <- data.frame(
  num_strains = integer(),
  num_private_orthogroups = integer()
)

for (i in 2:length(strainCol_c2_u)) {
  temp <- private_final %>% dplyr::select(2:i)
  # print(temp)
  
  private_count <- sum(apply(temp, 1, function(row) {if (sum(!is.na(row)) == 1) {sum(row, na.rm = TRUE)} else {0}}))
  rarefaction <- rbind(rarefaction, data.frame(num_strains = i-1, num_priv_genes = private_count))
}

rfc_genes <- ggplot(data = rarefaction) +
  geom_point(aes(x=num_strains, y = num_priv_genes), size=2, color = 'magenta3', alpha=0.8) +
  xlab("Strains") +
  ylab("Private genes") +
  # coord_cartesian(ylim = c(200,1500)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size=13, color = 'black', face = 'bold'),
    axis.text = element_text(size=11, color = 'black'),
    panel.border = element_rect(fill = NA),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  # 20pt on all sides
  )

rfc_genes

# COLOR points by Hawaiian strains or not!!!


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/rarefaction_115_labelled.png", height = 5, width = 9, rfc_genes, dpi = 600)


# Fit Michaelis-Menten model
mm_model <- nls(
  num_priv_genes ~ Vmax * num_strains / (Km + num_strains),
  data = rarefaction,
  start = list(Vmax = max(rarefaction$num_priv_genes), Km = median(rarefaction$num_strains))
)

# Extract parameters
Vmax_est <- coef(mm_model)["Vmax"]
Km_est <- coef(mm_model)["Km"]

# Define saturation proportion
p <- 0.95  # 95% saturation

# Analytical saturation estimate
x_sat <- (p / (1 - p)) * Km_est
threshold <- p * Vmax_est
additional_needed <- x_sat - max(rarefaction$num_strains)

# Predict over a reasonable x range for plotting
x_pred <- seq(min(rarefaction$num_strains), x_sat * 1.1, length.out = 1000)
y_pred <- predict(mm_model, newdata = data.frame(num_strains = x_pred))
pred_df <- data.frame(num_strains = x_pred, num_priv_genes = y_pred)

# Print summary
# cat("Estimated Vmax:", round(Vmax_est, 2), "\n")
# cat("Estimated Km:", round(Km_est, 2), "\n")
# cat("95% saturation reached at ~", round(x_sat, 2), "strains\n")
# cat("Additional strains needed:", round(additional_needed, 2), "\n")

rfc_genes_fit <- ggplot(rarefaction, aes(x = num_strains, y = num_priv_genes)) +
  geom_line(data = pred_df, aes(x = num_strains, y = num_priv_genes), color = "blue", size = 1) +
  geom_point(color = "magenta3", size = 2, alpha = 0.8) +
  geom_vline(xintercept = x_sat, linetype = "dashed", color = "black") +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
  annotate("text", x = x_sat * 0.95, y = threshold,
           label = paste0("Saturation at ~", round(x_sat), " strains"),
           hjust = 0.9, vjust = 3, color = "blue", size = 4) +
  xlab("Strains") +
  ylab("Private genes") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size=13, color = 'black', face = 'bold'),
    axis.text = element_text(size=11, color = 'black'),
    panel.border = element_rect(fill = NA),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  
  )
rfc_genes_fit

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/rarefaction_115_labelled_fit.png", height = 5, width = 9, rfc_genes_fit, dpi = 600)



# ======================================================================================================================================================================================== #
# Plotting distribution of strains and count of private genes #
# ======================================================================================================================================================================================== #














# ======================================================================================================================================================================================== #
# Plotting distribution of strains and count of private genes #
# ======================================================================================================================================================================================== #
trp <- private_final %>%
  dplyr::select(-Orthogroup) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("strain") %>%
  dplyr::mutate(strain = gsub("_count","", strain)) %>%
  dplyr::mutate(strain = gsub(".longest.protein_count","", strain)) %>%
  dplyr::mutate(count = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>%
  dplyr::select(strain, count)

trp$strain <- factor(trp$strain, levels = trp$strain)

ggplot(trp, aes(x = strain, y = count, fill = count)) +
  scale_fill_gradient(low = "pink", high = "magenta3") +
  geom_col() +
  theme(
    axis.text.x = element_text(color = 'black', angle = 70, hjust = 1),
    panel.background = element_blank(),
    axis.title = element_text(size=14, color = 'black'),
    axis.text.y = element_text(size=12, color = 'black'),
    panel.border = element_rect(fill = NA),
    legend.position = "none"
    ) +
  labs(x = "116 strains", y = "Private genes") 
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 90))

# Plot distribution of predicted PC genes for all strains
PC_count <- all_genes_strain %>%
  dplyr::count(strain, name = "n_genes") 

geneCount_HOGs <- trp %>%
  dplyr::left_join(PC_count, by = "strain")

geneCount_HOGs$strain <- factor(geneCount_HOGs$strain, levels = geneCount_HOGs$strain)
  
geneCount_HOGs <- geneCount_HOGs %>%
  dplyr::mutate(n_genes_scaled = n_genes / 100)

geneCount_HOGs$strain <- factor(geneCount_HOGs$strain, levels = geneCount_HOGs$strain)

HOGs_PCgenes <- ggplot(geneCount_HOGs, aes(x = strain)) +
  geom_col(aes(y = count), fill = "magenta3", width = 0.6) +
  geom_line(aes(y = n_genes_scaled, group = 1), color = "blue", size = 1.2) +
  geom_point(aes(y = n_genes_scaled), color = "blue", size = 2) +
  scale_y_continuous(
    expand = c(0.01, 0),
    name = "Private genes",
    sec.axis = sec_axis(~ ., name = "Protein-coding genes (1e02)")  
  ) +
  theme(
    axis.text.x = element_text(size = 6, color = 'black', angle = 70, hjust = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 22, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 18, color = 'black'),
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    # plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  # 20pt on all sides
  ) +
  labs(x = "Strains")
HOGs_PCgenes

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/HOGs_PCgenes.png", HOGs_PCgenes, height = 6, width = 10, dpi = 600)


lm_model <- lm(count ~ n_genes_scaled, data = geneCount_HOGs)
r2 <- summary(lm_model)$r.squared

blah <- ggplot(geneCount_HOGs, aes(x = n_genes, y = count)) +
  geom_point(color = 'magenta3') +
  geom_text(data = subset(geneCount_HOGs, n_genes > 27000 | n_genes < 20000 | count > 700), aes(label = strain), hjust = 0.9, vjust = 1.5, size = 5, color = "magenta3") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(r2, 3)), hjust = 2, vjust = 10, size = 7, color = "blue") +
  labs(y = "Private genes", x = "Protein-coding genes") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 22, color = 'black', face = 'bold'),
    axis.text = element_text(size = 18, color = 'black'),
    panel.border = element_rect(fill = NA))
blah

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/private_genes_PC.png", blah, height = 6, width = 10, dpi = 600)



# ======================================================================================================================================================================================== #
# Looking at relationship between number of private genes and number of single-exon genes #
# ======================================================================================================================================================================================== #

allGffs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/ALL_GFFs_longestIso_exons_fixed.tsv", col_names = c("exonID"))

singleEx <- allGffs %>%
  tidyr::separate(exonID, into = c("exonID","strain"), sep = " ") %>%
  dplyr::filter(strain != "c_elegans'") %>%
  dplyr::group_by(exonID, strain) %>%
  dplyr::mutate(num_exons = n()) %>%
  dplyr::filter(num_exons == 1) %>%
  dplyr::ungroup()

N2_exons <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/N2_longestIso_exons_fixed.tsv", col_names = c("exonID", "strain")) %>%
  dplyr::group_by(exonID,strain) %>%
  dplyr::mutate(num_exons = n()) %>%
  dplyr::filter(num_exons == 1) %>%
  dplyr::ungroup()

single_exonGenes <- singleEx %>%
  dplyr::bind_rows(N2_exons) %>%
  dplyr::count(strain, name = "num_SE_genes") %>%
  dplyr::filter(strain != "c_elegans")
  
priv_genes <- geneCount_HOGs %>%
  dplyr::select(strain,count, n_genes) %>%
  dplyr::rename(privates = count, PC_genes = n_genes) %>%
  dplyr::left_join(single_exonGenes, by = "strain") %>%
  dplyr::arrange(privates)


model <- lm(privates ~ num_SE_genes, data = priv_genes)
r_sq <- summary(model)$r.squared
segenes <- ggplot(priv_genes) + 
  geom_line(aes(x = privates, y = PC_genes/10), color = 'magenta') +
  geom_line(aes(x = privates, y = num_SE_genes), color = 'red') + 
  geom_point(aes(x = privates, y = num_SE_genes), color = 'firebrick') +
  geom_smooth(data = priv_genes, aes(x = privates, y = num_SE_genes),method = "lm", se = FALSE, color = "black") + 
  geom_text(data = subset(priv_genes, num_SE_genes > 3600 | num_SE_genes < 800), aes(x = privates, y = num_SE_genes, label = strain), hjust = 1.1, vjust = -0.5, size = 5, color = "firebrick") +
  annotate("text", x = 1000, y = 1500, label = "PC genes", hjust = -1.5, vjust = -5, size = 5, color = "magenta") +
  annotate("text", x = 1200, y = 4200, label = paste0("R² = ", round(r_sq, 3)), size = 5, color = "black") +
  xlab("Number of private genes") +
  ylab("Single exon genes") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size=13, color = 'black', face = 'bold'),
    axis.text = element_text(size=11, color = 'black'),
    panel.border = element_rect(fill = NA),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  
  )
segenes


# Of the private genes, are all of them single exon genes??? 
private_ogs <- private %>%
  dplyr::select(HOG) %>%
  dplyr::pull()


N2 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/N2_tran_gene.tsv", col_names = c('tran','gene'))

strain_SE <- singleEx %>%
  dplyr::bind_rows(N2_exons) %>%
  dplyr::select(-num_exons) %>%
  dplyr::filter(strain != "c_elegans") %>%
  dplyr::mutate(exonID = gsub("\\.t[0-9]+", "", exonID)) %>%
  dplyr::mutate(exonID = gsub("transcript:", "", exonID)) %>%
  dplyr::mutate(exonID = ifelse(exonID %in% N2$tran, N2$gene, exonID)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_single_exon = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(strain, exonID, num_single_exon)
  
privs <- ortho_genes_dd %>% 
  dplyr::select(-OG,-'Gene Tree Parent Clade') %>%
  dplyr::bind_rows((private_OGs %>% dplyr::rename(HOG = Orthogroup))) %>%
  dplyr::filter(HOG %in% private_ogs) %>%
  pivot_longer(
    cols = -HOG,
    names_to = "strain",
    values_to = "genes"
  ) %>%
  dplyr::filter(!is.na(genes)) %>%                     
  dplyr::mutate(genes = str_split(genes, ",")) %>%     
  unnest(genes) %>%                             
  dplyr::arrange(strain) %>%
  dplyr::select(-HOG) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_privs = n()) %>%
  dplyr::ungroup() 

prop_SE_priv <- privs %>%
  inner_join(strain_SE, by = c("strain" = "strain", "genes" = "exonID")) %>%  
  group_by(strain) %>%
  summarise(
    num_priv_single_exon = n(),                        
    num_privs = first(num_privs),                      
    proportion = num_priv_single_exon / num_privs      
  )

plot_data <- prop_SE_priv %>%
  dplyr::mutate(non_single_exon = 1 - proportion) %>%
  dplyr::arrange(proportion) %>% 
  dplyr::mutate(strain = factor(strain, levels = strain)) %>%
  dplyr::select(strain, single_exon = proportion, non_single_exon) %>%
  pivot_longer(cols = c(single_exon, non_single_exon), 
               names_to = "gene_type", values_to = "prop") 

SE_plot <- ggplot(plot_data, aes(x = prop, y = reorder(strain, prop), fill = gene_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("single_exon" = "firebrick", "non_single_exon" = "gray70"),
                    labels = c("multi-exon", "single-exon")) +
  labs(x = "Proportion of private genes", y = "Strain") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9, color = 'black'),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.border = element_rect(fill = NA)) +
  scale_x_continuous(expand = c(0,0))
SE_plot




# SEE SCRIPT "pangenome_Ce_geneSet_GO.R
# # ======================================================================================================================================================================================== #
# # Gene Ontology of N2 genes in each gene set #
# # ======================================================================================================================================================================================== #
# # Extract HOGs that have ONE N2 gene to lift over ontology to N2 orthologs #
# # all_relations_classification_rowid <- count %>% dplyr::rename(rowid = HOG)
# # all_relations_rowid <- all_relations %>% dplyr::rename(rowid = HOG)
# # ortho_genes_dd_rowid <- ortho_genes_dd %>% dplyr::rename(rowid = HOG)
# 
# # Extract rowids for HOGs that have one N2 gene contributing
# oneN2 <- all_relations %>%
#   dplyr::filter(N2_count == 1) %>% 
#   dplyr::pull(HOG)
# 
# private_n2 <- count %>%
#   dplyr::filter(class == "private") %>%
#   dplyr::filter(!is.na(N2_count)) %>%
#   dplyr::pull(HOG) # 16 - this is correct
# 
# one_n2_and_private <- c(oneN2,private_n2)
# 
# HOG_classification <- count %>%
#   dplyr::select(HOG,freq,class)
# 
# n2_genes_GO_liftover <- ortho_genes_dd %>% dplyr::filter(HOG %in% one_n2_and_private) %>%
#   dplyr::select(-OG, -"Gene Tree Parent Clade") %>%
#   dplyr::left_join(HOG_classification, by = "HOG") %>%
#   tidyr::separate_rows(N2, sep = ",\\s*") # need to split N2 genes in private HOGs because there are multiple N2 genes 
# 
# all_N2_genes <- all_genes_strain %>%
#   dplyr::filter(strain == "N2") %>%
#   dplyr::pull(attributes)
# 
# genes_core <- n2_genes_GO_liftover %>%
#   dplyr::filter(class == "core") %>%
#   dplyr::pull(N2) # 14,550
# 
# genes_acc <- n2_genes_GO_liftover %>%
#   dplyr::filter(class == "accessory") %>%
#   dplyr::pull(N2) # 4,236
#   
# genes_private <- n2_genes_GO_liftover %>%
#   dplyr::filter(class == "private") %>%
#   dplyr::pull(N2) # 37
# 
# 
# enrich_go <- function(wb_ids){
#   
#   mf <- enrichGO(gene          = wb_ids,
#                  OrgDb         = org.Ce.eg.db,
#                  ont           = "MF",
#                  pAdjustMethod = "bonferroni",
#                  keyType       = 'ENSEMBL',
#                  pvalueCutoff  = 0.05,
#                  qvalueCutoff  = 0.05)
#   
#   bp <- enrichGO(gene          = wb_ids,
#                  OrgDb         = org.Ce.eg.db,
#                  ont           = "BP",
#                  pAdjustMethod = "bonferroni",
#                  keyType       = 'ENSEMBL',
#                  pvalueCutoff  = 0.05,
#                  qvalueCutoff  = 0.05)
#   
#   return(list(mf, bp))
# }
# 
# 
# 
# N2_anno <- enrich_go(all_N2_genes)
# core_anno <- enrich_go(genes_core)
# acc_anno <- enrich_go(genes_acc)
# private_anno <- enrich_go(genes_private)
# 
# data_man_plot <- function(N2_anno, geneSet_anno, gene_set) {
#   genesBP_N2 <- setReadable(N2_anno[[2]], OrgDb = org.Ce.eg.db)
#   gene_BP_geneSet <- setReadable(geneSet_anno[[2]], OrgDb = org.Ce.eg.db)
#   
#   df_GO_enrich_BP <- rbind(dplyr::mutate(genesBP_N2[], freq="control_N2_genes"), dplyr::mutate(gene_BP_geneSet[], freq = paste0(gene_set,"_genes"))) %>%
#     tidyr::separate(GeneRatio, c("genes_enrich", "genes_in_database")) %>%
#     tidyr::separate(BgRatio, c("genes_in_geneSet", "N2_genes_in_database")) %>%
#     dplyr::mutate(genes_enrich = as.numeric(genes_enrich), genes_in_database = as.numeric(genes_in_database), genes_in_geneSet = as.numeric(genes_in_geneSet), N2_genes_in_database = as.numeric(N2_genes_in_database)) %>%
#     dplyr::mutate(GeneRatio = genes_enrich/genes_in_database, BgRatio = genes_in_geneSet/N2_genes_in_database) %>%
#     dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
#     dplyr::arrange(p.adjust)
#   
#   GO_list_BP <- df_GO_enrich_BP$Description
#   
#   GO_list_BP_plotpoint <- data.frame(Description=GO_list_BP, plotpoint=length(GO_list_BP):1)
#   
#   df_GO_enrich_BP_sum <- df_GO_enrich_BP %>%
#     dplyr::filter(Description %in% GO_list_BP) %>%
#     dplyr::left_join(., GO_list_BP_plotpoint, by = "Description") %>%
#     dplyr::group_by(ID) %>%
#     dplyr::mutate(class_gene_total = max(genes_in_geneSet)) %>%
#     dplyr::ungroup()
#   
#   df_class_total_BP <- df_GO_enrich_BP_sum %>%
#     dplyr::distinct(Description, class_gene_total, plotpoint) %>%
#     dplyr::arrange(-plotpoint)
#   
#   df_GO_enrich_BP_sum$freq <- factor(df_GO_enrich_BP_sum$freq, levels = c(paste0(gene_set,"_genes"),"control_N2_genes"), labels = c(gene_set,"All N2 genes"))
#   
#   if (gene_set == "Core") { 
#     plot_GO_BP <- ggplot(df_GO_enrich_BP_sum) +
#       geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
#       geom_point(aes(x = -log10(p.adjust), y = plotpoint, shape = freq, size = EnrichRatio, fill = Count)) +
#       scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.02,0.02)) +
#       scale_fill_gradient(low = "darkolivegreen1", high = "green4", breaks = c(round(min(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(median(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(max(df_GO_enrich_BP_sum$Count, na.rm = TRUE)))) +
#       scale_shape_manual(values=c(21, 22)) +
#       scale_size_continuous(range = c(1,4)) + #breaks = c(1, 2, 3)) + 
#       theme(axis.text.x = element_text(size=9, color='black'), 
#             axis.text.y = element_text(size=7, color='black'),
#             axis.title = element_text(size=10, color='black', face = 'bold'), 
#             legend.title = element_text(size=7, color='black'), 
#             legend.text = element_text(size=6, color='black'), 
#             legend.position = "inside", 
#             legend.position.inside = c(0.75, 0.18),
#             plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#             panel.grid = element_blank(),
#             panel.background = element_blank(),
#             panel.border = element_rect(fill = NA),
#             legend.direction = "horizontal", legend.box = "vertical",
#             plot.margin = margin(l = 20, unit = "pt"),
#             legend.spacing.y = unit(0.02, 'in'),
#             # legend.spacing.x = unit(0.02, "cm"),
#             legend.key.size = unit(0.3, "cm")) +
#       guides(shape = guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3), force = TRUE), 
#              size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE), 
#              fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2, force = TRUE)) +
#       labs(title = paste0(gene_set," gene set"), x="-log10 (corrected p-value)", shape = "Gene Set", size = "Fold enrichment", fill = "Gene counts") 
#   } else if (gene_set == "Accessory") {
#     plot_GO_BP <- ggplot(df_GO_enrich_BP_sum) +
#       geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
#       geom_point(aes(x = -log10(p.adjust), y = plotpoint, shape = freq, size = EnrichRatio, fill = Count)) +
#       scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.02,0.02)) +
#       scale_fill_gradient(low = "burlywood1", high = "#DB6333", breaks = c(round(min(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(median(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(max(df_GO_enrich_BP_sum$Count, na.rm = TRUE)))) +
#       scale_shape_manual(values=c(21, 22)) +
#       scale_size_continuous(range = c(2,10), breaks = c(3, 3.75, 4.5)) + 
#       theme(axis.text.x = element_text(size=9, color='black'), 
#             axis.text.y = element_text(size=7, color='black'),
#             axis.title = element_text(size=10, color='black', face = 'bold'), 
#             legend.title = element_text(size=7, color='black'), 
#             legend.text = element_text(size=6, color='black'), 
#             plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#   
#             legend.position = "inside", 
#             legend.position.inside = c(0.72, 0.18),
#             panel.grid = element_blank(),
#             panel.background = element_blank(),
#             panel.border = element_rect(fill = NA),
#             legend.direction = "horizontal", legend.box = "vertical",
#             # plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
#             legend.spacing.y = unit(0.02, 'in'),
#             # legend.spacing.x = unit(0.02, "cm"),
#             legend.key.size = unit(0.3, 'cm')) +
#       guides(shape = guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3), force = TRUE), 
#              size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE), 
#              fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2, force = TRUE)) +
#       labs(title = paste0(gene_set," gene set"), x="-log10 (corrected p-value)", shape = "Gene Set", size = "Fold enrichment", fill = "Gene counts") 
#   } else if (gene_set == "Private") {
#     plot_GO_BP <- ggplot(df_GO_enrich_BP_sum) +
#       geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
#       geom_point(aes(x = -log10(p.adjust), y = plotpoint, shape = freq, size = EnrichRatio, fill = Count)) +
#       scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.02,0.02)) +
#       scale_fill_gradient(low = "thistle1", high = "magenta3", breaks = c(round(min(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(median(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(max(df_GO_enrich_BP_sum$Count, na.rm = TRUE)))) +
#       scale_shape_manual(values=c(21, 22)) +
#       scale_size_continuous(range = c(2,8), breaks = c(50, 150, 250)) + 
#       theme(axis.text.x = element_text(size=9, color='black'), 
#             axis.text.y = element_text(size=7, color='black'),
#             axis.title = element_text(size=10, color='black', face = 'bold'), 
#             legend.title = element_text(size=7, color='black'), 
#             legend.text = element_text(size=6, color='black'), 
#             legend.position = "inside", 
#             legend.position.inside = c(0.72, 0.18),
#             panel.grid = element_blank(),
#             plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
#             panel.background = element_blank(),
#             panel.border = element_rect(fill = NA),
#             legend.direction = "horizontal", legend.box = "vertical",
#             # plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
#             legend.spacing.y = unit(0.02, 'in'),
#             #legend.spacing.x = unit(0.02, "cm"),
#             legend.key.size = unit(0.3, 'cm')) +
#       guides(shape = guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3), force = TRUE), 
#              size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE), 
#              fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2, force = TRUE)) +
#       labs(title = paste0(gene_set," gene set"), x="-log10 (corrected p-value)", shape = "Gene Set", size = "Fold enrichment", fill = "Gene counts") 
#   } else {
#     stop("Please provide gene set as 'Core', 'Accessory', or 'Private'")
#   }
#   
#   return(plot_GO_BP)
# }
# 
# core_plot <- data_man_plot(N2_anno, core_anno, "Core")
# core_plot
# 
# acc_plot <- data_man_plot(N2_anno, acc_anno, "Accessory")
# acc_plot
# 
# private_plot <- data_man_plot(N2_anno, private_anno, "Private")
# private_plot
# 
# 
# 
# GO_all <- plot_grid(
#   core_plot, acc_plot, private_plot,
#   nrow = 1,
#   rel_heights = c(1, 1), rel_widths = c(1,1,1),
#   align = "hv"
# )
# GO_all
# 
# 
# # ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/GO_pangenomeGeneSet.png", GO_all, height = 8, width = 26, dpi = 600)
