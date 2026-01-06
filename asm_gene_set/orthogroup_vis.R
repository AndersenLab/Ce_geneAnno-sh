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
library(org.Ce.eg.db)
library(GO.db)
library(clusterProfiler) ## BiocManager::install("clusterProfiler") # need this and the next package???
library(enrichplot)


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
# genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
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

  # dplyr::mutate(gene = str_extract(attributes,"(?<=\\bID=gene:)WBGene\\d+|^WBGene\\d+|(?<=\\bName=)WBGene\\d+"),tran = str_extract(attributes, "(?<=\\bsequence_name=)[^;]+")) %>%
  # dplyr::select(tran,gene)

# test1 <- all_genes_strain %>%
#   dplyr::filter(strain == "N2") # ~19,000

# write.table(all_genes_strain,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/116strain_genes.tsv", quote = F, row.names = F, col.names = T, sep = '\t')

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-IDY,-LENR,-LENQ) %>% dplyr::filter(strain != "ECA396")
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
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))
# whatever <- ortho_genes_dd %>%
  # dplyr::select(JU1581.braker.longest.protein, N2.longest.protein)

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

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/gene_set_allOrtho_142.png", gs_allOrtho, height = 5, width = 11, dpi = 600)


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



# ======================================================================================================================================================================================== #
# PLOTTING HORIZONTAL BAR PLOTS FOR THE PROPORTION OF GENES THAT ARE CLASSIFIED AS EACH GENE SET IN STRAIN
# ======================================================================================================================================================================================== #
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

N2_gene_count <- df_long %>% dplyr::filter(strain == "N2")

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
  geom_bar(stat = "identity", alpha = 0.5, color = "black", linewidth = 0.3) +
  # geom_vline(xintercept = core$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
  # geom_vline(xintercept = acc$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
  # geom_vline(xintercept = priv$meann, linetype = "dashed", color = 'gray22', linewidth = 1.5) +
  labs(x = "Percent of genes", fill = "Gene set") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    axis.text.y = element_text(size = 7.5, color = 'black'),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 22, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)
  )
contrib

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/geneSet_contrib_142.png", contrib, height = 13, width = 12, dpi = 600)







# ======================================================================================================================================================================================== #
# Pangenome gene set classification of N2 genes enriched in particular N2 genomic loci??
# ======================================================================================================================================================================================== #
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

N2_coords <- N2_gff %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(attributes = gsub("ID=transcript:","",attributes), attributes = gsub("Parent=gene:","",attributes)) %>%
  tidyr::separate_wider_delim(attributes, delim = ";",names = c("tran", "gene", "rest"), too_many = "merge") %>%
  dplyr::select(seqid, start, end, tran, gene, -rest)

n2_gene <- N2_coords %>%
  dplyr::mutate(tran = paste0("transcript_",tran))

n2_table <- ortho_genes_dd %>%
  dplyr::bind_rows(private_OGs) %>% 
  dplyr::select(Orthogroup,N2)

ortho_count_wCoord <- count %>%
  dplyr::left_join(n2_table, by = "Orthogroup") %>%
  dplyr::select(freq, class, N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ",\\s*") %>% # splitting rows so each gene is on a row and it retains is gene set classification
  dplyr::left_join(n2_gene, by = c("N2" = "tran")) %>%
  dplyr::select(freq,class,N2,gene,seqid,start,end)
  

plt_data <- ortho_count_wCoord %>%
  dplyr::filter(seqid != "MtDNA") %>%
  dplyr::mutate(mid_mb = (start + end) / 2 / 1e6)
  
ggplot(data = plt_data) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, fill = class), ymin = -Inf, ymax = Inf, alpha = 0.5) +
  geom_density(aes(x = mid_mb, y = after_stat(scaled), color = class), adjust = 0.5, linewidth = 0.9, position = "identity") +
  scale_color_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3", "core" = "green4")) +
  facet_wrap(~seqid, scales = "free") +
  theme(axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
        axis.text.x = element_text(size = 13, color = 'black'),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("N2 genome position (Mb)") +
  scale_y_continuous(NULL, breaks = NULL) 


plt_data <- plt_data %>%
  mutate(seqid = factor(seqid, levels = c("I","II","III","IV","V","X"))) %>%
  dplyr::rename(Class = class)
ncol_facets <- 3
lev <- levels(plt_data$seqid)
left_facets <- lev[seq(1, length(lev), by = ncol_facets)] 

class_labels <- plt_data %>% dplyr::distinct(Class) %>% dplyr::mutate(y = ifelse(Class == "core", 1.75,
                                                                                 ifelse(Class == "accessory", 1, 0.25)), x = -Inf)
class_labels <- tidyr::crossing(class_labels, seqid = left_facets)

ggplot() + 
  geom_rect(data = plt_data %>% dplyr::filter(Class == "core"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 1.5, ymax = 2), fill = "green4", alpha = 0.5) +
  geom_rect(data = plt_data %>% dplyr::filter(Class == "accessory"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.75, ymax = 1.25), fill = "#DB6333", alpha = 0.5) +
  geom_rect(data = plt_data %>% dplyr::filter(Class == "private"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 0.5), fill = "magenta3", alpha = 0.5) +
  geom_text(data = class_labels,aes(x = x, y = y, label = Class), hjust = 1.1, fontface = "bold", size = 4, inherit.aes = FALSE) +
  facet_wrap(~seqid, scales = "free") +
  coord_cartesian(clip = "off") +
  theme(axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
        axis.text.x = element_text(size = 13, color = 'black'),
        panel.background = element_blank(),
        panel.border = element_rect(fill = "NA", color = 'black'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(l = 70)) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("N2 genome position (Mb)") +
  scale_y_continuous(NULL, breaks = NULL)


ggplot(data = plt_data) +
  geom_density(aes(x = mid_mb, y = after_stat(count), color = Class), adjust = 1, linewidth = 0.9, position = "identity") + 
  scale_color_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  facet_wrap(~seqid, scales = 'free_x') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = "NA", color = 'black'),
        # strip.background = element_blank(),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title.y = element_text(size =14, color = 'black', face = 'bold'),
        # strip.text.x = element_blank(),
        plot.margin = margin(l = 70)) +
  ylab("Kernel Density * Count")


# Final plot
lane_h   <- 14        # lane height (y units)
gap_h    <- 6         # gap between lanes
y_core   <- -lane_h                      # [-12, 0)
y_acc    <- -(2*lane_h + gap_h)          # [-30, -18)
y_priv   <- -(3*lane_h + 2*gap_h)        # [-48, -36)

rects <- plt_data %>%
  dplyr::mutate(xmin = start/1e6, xmax = end/1e6, ymin = dplyr::case_when(
    Class == "core" ~ y_core - lane_h, 
    Class == "accessory" ~ y_acc - lane_h, 
    Class == "private" ~ y_priv - lane_h), 
    ymax = dplyr::case_when(
      Class == "core" ~ y_core, 
      Class == "accessory" ~ y_acc, 
      Class == "private" ~ y_priv)) %>% 
  dplyr::rename(`Gene set` = Class)

class_labels_final <- plt_data %>% dplyr::distinct(Class) %>% dplyr::mutate(y = ifelse(Class == "core", -21,
                                                                                       ifelse(Class == "accessory", -41, -61)), x = -Inf)
class_labels_final <- tidyr::crossing(class_labels_final, seqid = left_facets)

density_label <- plt_data %>% dplyr::mutate(label = "Gene count density") %>% dplyr::mutate(x = -Inf, y = 100) %>% dplyr::select(label,x,y)
density_label_final <- tidyr::crossing(density_label, seqid = left_facets)

finalfinal <- ggplot() +
  # draw rect "tracks" first so they sit under the density
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = `Gene set`), color = NA) +
  scale_fill_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  # density scaled by counts (absolute abundance)
  geom_density(data = plt_data, aes(x = mid_mb, y = after_stat(count), color = Class), adjust = 0.5, linewidth = 0.9, position = "identity", show.legend = FALSE) +
  geom_text(data = class_labels_final, aes(x = x, y = y, label = Class), hjust = 1.1, fontface = "bold", size = 4, inherit.aes = FALSE) +
  geom_text(data = density_label_final, aes(x = x, y = y, label = label), vjust = -3, fontface = "bold", size = 5, inherit.aes = FALSE, angle = 90) +
  scale_color_manual(values = c(accessory="#DB6333", private="magenta3", core="green4")) +
  facet_wrap(~ seqid, scales = "free_x", ncol = 3) +
  # reserve room below zero for the lanes
  coord_cartesian(ylim = c(y_priv - 5, NA), clip = "off") +
  scale_x_continuous(expand = c(0.01,0)) +
  # hide negative tick labels; show only non-negative y ticks
  scale_y_continuous(name = "Kernel Density × Count", labels = function(b) ifelse(b < 0, "", b)) + # normalized kernel density * count
  labs(x = "N2 genome position (Mb)") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    # strip.background = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text = element_text(size =10, color = 'black'),
    axis.text.x  = element_text(size = 12),
    panel.spacing=unit(2, "lines"),
    # axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = 'none',
    axis.text.y  = element_text(size = 12),
    plot.margin  = margin(l = 70, r = 5, t = 5, b = 5))
finalfinal

bw_by_class <- plt_data %>%
  dplyr::group_by(Class) %>%
  dplyr::summarise(bw = density(mid_mb)$bw)
bw_by_class


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/N2_geneSetLoci_density.png",finalfinal, height = 12, width = 14, dpi = 500)

# 
# 
#   
# ortho_private_coords <- ortho_count_wCoord %>% dplyr::filter(class == "private")
# 
# summ <- ortho_private_coords %>%
#   dplyr::group_by(seqid) %>%
#   dplyr::summarise(n_perChrom = n())
# 
# ggplot(data = ortho_private_coords) +
#   geom_rect(data = N2_coords, aes(xmin = start / 1e6, xmax = end / 1e6), ymin = -Inf, ymax = Inf, fill = 'gray30', alpha = 0.5) + 
#   geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6), ymin = -Inf, ymax = Inf, fill = 'magenta3', alpha = 0.5) + 
#   # scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3", "core" = "green4")) +
#   facet_wrap(~seqid, scales = "free") +
#   theme(axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
#         axis.text.x = element_text(size = 13, color = 'black'),
#         panel.background = element_blank(),
#         panel.border = element_rect(fill = NA, color = 'black')) +
#   scale_x_continuous(expand = c(0,0)) +
#   xlab("N2 genome position (Mb)")
# 



# ======================================================================================================================================================================================== #
# PLOTTING PANGENOME AND CORE PANGENOME RAREFACTION CURVES # 
# ======================================================================================================================================================================================== #
# Randomly sample 1-141 strains and calculate the number of total orthogroups, and the number of core orthogroups
all <- all_relations %>% dplyr::select(-Orthogroup)

pan_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/pan_iterativeOGcount.rds")
core_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/core_iterativeOGcount.rds")
accessory_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/acc_iterativeOGcount.rds")
priv_final <- readRDS("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/priv_iterativeOGcount.rds")

set.seed(42)

n_strains_total <- ncol(all)
n_perms <- 100
 
 
# # For the pangenome
# pan_list <- vector("list", length = n_strains_total - 1) # -1 because we iterate 2 - 142
# iteration_pan <- 1
#
# for (i in 2:n_strains_total) {
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     subset <- subset %>% dplyr::filter(!if_all(everything(), is.na))
#     all_OGs <- nrow(subset)
#     # print(all_OGs)
#     print(paste0("On strain subset: ",i,", and iteration: ", it_i))
# 
#     pan_list[[iteration_pan]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_core_ogs = all_OGs)
#     iteration_pan <- iteration_pan + 1
#   }
# }
# 
# pan_final <- dplyr::bind_rows(pan_list)

pan_summary <- pan_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_core = median(n_core_ogs),
    mean_core   = mean(n_core_ogs),
    sd_core     = sd(n_core_ogs),
    q05         = quantile(n_core_ogs, 0.05),
    q95         = quantile(n_core_ogs, 0.95)) %>%
  dplyr::ungroup()

# # For the core pangenome
# res_list <- vector("list", length = n_strains_total - 1) # -1 because we iterate 2 - 141
# iteration <- 1
# 
# for (i in 2:n_strains_total) {
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     # binarize and count core
#     core_calc <- subset %>%
#       dplyr::mutate(across(everything(), ~ ifelse(is.na(.),0, ifelse(. >= 1, 1, .)))) %>%
#       dplyr::mutate(sum = rowSums(across(everything()))) %>%
#       dplyr::mutate(freq = (sum / i)) %>%
#       dplyr::mutate(class = case_when(freq == 1 ~ "core")) %>%
#       dplyr::filter(class == "core")
# 
#     # print(head(core_calc))
#     print(paste0("On strain subset: ", i,", and iteration: ", it_i))
# 
#     core_count <- nrow(core_calc)
#     # print(core_count)
# 
#     res_list[[iteration]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_core_ogs = core_count)
#     iteration <- iteration + 1
#   }
# }
# 
# core_final <- dplyr::bind_rows(res_list)

core_summary <- core_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_core = median(n_core_ogs),
    mean_core   = mean(n_core_ogs),
    sd_core     = sd(n_core_ogs),
    q05         = quantile(n_core_ogs, 0.05),
    q95         = quantile(n_core_ogs, 0.95)) %>%
  dplyr::ungroup()


# For the accessory pangenome
# res_list <- vector("list", length = n_strains_total - 1) # -1 because we iterate 2 - 141
# iteration <- 1
# 
# for (i in 2:n_strains_total) {
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     # binarize and count accessory
#     accessory_calc <- subset %>%
#       dplyr::mutate(across(everything(), ~ ifelse(is.na(.),0, ifelse(. >= 1, 1, .)))) %>%
#       dplyr::mutate(sum = rowSums(across(everything()))) %>%
#       dplyr::mutate(freq = (sum / i)) %>%
#       dplyr::mutate(
#         class = case_when(
#           freq == 1 ~ "core",
#           freq > 1/i & freq < 1 ~ "accessory",
#           freq == 1/i ~ "private",
#           TRUE ~ "undefined")) %>%
#       dplyr::filter(class == "accessory")
# 
#     # print(head(accessory_calc))
#     print(paste0("On strain subset: ", i,", and iteration: ", it_i))
# 
#     accessory_count <- nrow(accessory_calc)
#     # print(accessory_count)
# 
#     res_list[[iteration]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_accessory_ogs = accessory_count)
#     iteration <- iteration + 1
#   }
# }
# 
# accessory_final <- dplyr::bind_rows(res_list)

accessory_summary <- accessory_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_accessory = median(n_accessory_ogs),
    mean_accessory   = mean(n_accessory_ogs),
    sd_accessory     = sd(n_accessory_ogs),
    q05         = quantile(n_accessory_ogs, 0.05),
    q95         = quantile(n_accessory_ogs, 0.95)) %>%
  dplyr::ungroup()

# For the private pangenome
# res_list <- vector("list", length = n_strains_total -1)
# iteration <- 1
# 
# for (i in 2:n_strains_total) {              # NEED TO BEGIN ITERATION ON A SINLE STRAIN WHEN CALCULATING PRIVATE
#   for (it_i in 1:n_perms) {
#     # pick k random strains
#     cols <- sample(colnames(all), size = i, replace = FALSE)
#     subset <- all[, cols, drop = FALSE] # subset all df to k strains
# 
#     # binarize and count accessory
#     priv_calc <- subset %>%
#       dplyr::mutate(across(everything(), ~ ifelse(is.na(.),0, ifelse(. >= 1, 1, .)))) %>%
#       dplyr::mutate(sum = rowSums(across(everything()))) %>%
#       dplyr::mutate(freq = (sum / i)) %>%
#       dplyr::mutate(
#         class = case_when(
#           freq == 1 ~ "core",
#           freq > 1/i & freq < 1 ~ "accessory",
#           freq == 1/i ~ "private",
#           TRUE ~ "undefined")) %>%
#       dplyr::filter(class == "private")
# 
#     # print(head(priv_calc))
#     print(paste0("On strain subset: ", i,", and iteration: ", it_i))
# 
#     priv_count <- nrow(priv_calc)
#     # print(priv_count)
# 
#     res_list[[iteration]] <- data.frame(
#       n_strains = i,
#       replicate = it_i,
#       n_priv_ogs = priv_count)
#     iteration <- iteration + 1
#   }
# }
# 
# priv_final <- dplyr::bind_rows(res_list)

priv_summary <- priv_final %>%
  dplyr::group_by(n_strains) %>%
  dplyr::summarise(
    median_priv = median(n_priv_ogs),
    mean_priv   = mean(n_priv_ogs),
    sd_priv     = sd(n_priv_ogs),
    q05         = quantile(n_priv_ogs, 0.05),
    q95         = quantile(n_priv_ogs, 0.95)) %>%
  dplyr::ungroup()



# Plotting
pan_rarefact <- ggplot() +
  # Pangenome
  geom_errorbar(data = pan_summary, aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), width = 0.5) +
  geom_point(data = pan_summary, aes(x = n_strains, y = mean_core), color = 'blue', size = 3) +
  # Core
  geom_errorbar(data = core_summary, aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), width = 0.5) +
  geom_point(data = core_summary, aes(x = n_strains, y = mean_core), color = 'green4', size = 3) +
  # Accessory
  geom_errorbar(data = accessory_summary, aes(x = n_strains, ymin = mean_accessory - sd_accessory, ymax = mean_accessory + sd_accessory), width = 0.5) +
  geom_point(data = accessory_summary, aes(x = n_strains, y = mean_accessory), color = '#DB6333', size = 3) +
  # Private
  geom_errorbar(data = priv_summary, aes(x = n_strains, ymin = mean_priv - sd_priv, ymax = mean_priv + sd_priv), width = 0.5) +
  geom_point(data = priv_summary, aes(x = n_strains, y = mean_priv), color = 'magenta3', size = 3) +
  # geom_ribbon(aes(x = n_strains, ymin = mean_core - sd_core, ymax = mean_core + sd_core), alpha = 0.2) +
  labs(x = "Genomes", y = "Orthogroups") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size =14, color = 'black')) +
  coord_cartesian(ylim = c(0,60000))
    # legend.position = 'none')
pan_rarefact

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/pan_core_rarefaction.png", pan_rarefact, width = 14, height = 12, dpi = 600)


# 
# saveRDS(pan_final, file = "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/pan_iterativeOGcount.rds")
# saveRDS(core_final, file = "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/core_iterativeOGcount.rds")
# saveRDS(accessory_final, file = "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/acc_iterativeOGcount.rds")
# saveRDS(priv_final, file = "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/plots/priv_iterativeOGcount.rds")









# ======================================================================================================================================================================================== #
# Plotting orthogroup/gene family presence/absence heatmap, clustered by strain relatedness
# ======================================================================================================================================================================================== #
colnames(all_relations) <- strainCol_c2
pav_mat <- all_relations %>% dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>% 
  dplyr::mutate(across(2:ncol(.), ~ifelse(is.na(.), 0, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(-freq) %>%
  dplyr::mutate(class = factor(class, levels = c("core", "accessory", "private"))) %>%
  dplyr::arrange(desc(sum)) %>%
  dplyr::select(-sum) %>%
  dplyr::mutate(Orthogroup = factor(Orthogroup, levels = unique(Orthogroup)))

pav_long <- pav_mat %>%
  tidyr::pivot_longer(
    cols = all_of(strainCol_c2_u),
    names_to = "strain",
    values_to = "presence") %>%
  dplyr::mutate(
    class_presence = case_when(
      presence == 1 ~ paste0(class, "_present"),
      presence == 0 ~ paste0(class, "_absent"))) %>%
  dplyr::group_by(strain,class) %>%
  dplyr::mutate(class_count=sum(presence)) %>%
  dplyr::ungroup() 

strain_order_acc <- pav_long %>%
  dplyr::filter(class=="accessory") %>%
  dplyr::arrange(class_count) %>%
  dplyr::distinct(strain,.keep_all = T) %>%
  dplyr::pull(strain) 

ggplot(pav_long %>%
         dplyr::mutate(strain=factor(strain,levels=strain_order_acc)),
       aes(x = Orthogroup, y = strain, fill = class_presence)) +
  geom_tile() +
  scale_fill_manual(
    values = c(
      core_present       = "green4",
      core_absent        = "white",
      accessory_present  = "#DB6333",
      accessory_absent   = "white",
      private_present    = "magenta3",
      private_absent     = "white")) +
  labs(x = "Orthogroups", y = "Strain", fill = "Class × presence") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),      
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 6, color = 'black'),
    axis.title.x = element_text(color = 'black', size = 12, face = 'bold'),
    axis.title.y = element_blank())














# ======================================================================================================================================================================================== #
# Plotting N2 genes frequency with HDR resolution
# ======================================================================================================================================================================================== #
# Turning gene count table into a binary presence/absense of an orthogroup 
# count <- all_relations %>%
#   dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
#   dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
#   dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
#   dplyr::mutate(
#     class = case_when(
#       freq == 1 ~ "core",
#       freq > private_freq & freq < 1 ~ "accessory",
#       freq == private_freq ~ "private",
#       TRUE ~ "undefined"
#     )
#   ) 
# 
# 
# n2_gene <- all_genes_strain %>%
#   dplyr::filter(strain == "N2") %>%
#   dplyr::rename(ID = attributes)
# 
# n2_table <- ortho_genes_dd %>% 
#   dplyr::select(-OG,-'Gene Tree Parent Clade') %>%
#   dplyr::bind_rows((private_OGs %>% dplyr::rename(HOG = Orthogroup))) %>%
#   dplyr::select(HOG,N2) 
# 
# ortho_count_wCoord <- count %>%
#   dplyr::left_join(n2_table, by = "HOG") %>%
#   dplyr::select(freq, class, N2) %>%
#   dplyr::filter(!is.na(N2)) %>%
#   tidyr::separate_rows(N2, sep = ",\\s*") %>% # splitting rows so each gene is on a row and it retains is gene set classification
#   dplyr::left_join(n2_gene, by = c("N2" = "ID")) %>%
#   dplyr::select(freq,class,N2,contig,start,end) 
# 
# HDRs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/114HDRs.tsv", col_names = c("CHROM","START","END","strain"))
# 
# all_regions <- HDRs %>%
#   dplyr::arrange(CHROM,START) %>%
#   dplyr::group_split(CHROM)
# 
# strain_count <- HDRs %>% dplyr::distinct(strain, .keep_all = T)
# print(nrow(strain_count))
# 
# 
# getRegFreq <- function(all_regions) {
#   all_collapsed <- list()
#   for (i in 1:length(all_regions)) {
#     temp <- all_regions[[i]]
#     k=1
#     j=1
#     while (k==1) {
#       print(paste0("chrom:",i,"/iteration:",j))
#       checkIntersect <- temp %>% 
#         dplyr::arrange(CHROM,START) %>%
#         dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
#         dplyr::mutate(check=ifelse(is.na(check),F,check))
#       
#       #print(nrow(checkIntersect %>% dplyr::filter(check==T)))
#       
#       if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
#         print("NO MORE INTERSECTS")
#         k=0
#       } else {
#         
#         temp <- checkIntersect %>%
#           dplyr::mutate(gid=data.table::rleid(check)) %>%
#           dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))
#         
#         collapse <- temp %>%
#           dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
#           dplyr::group_by(gid) %>%
#           dplyr::mutate(newStart=min(START)) %>%
#           dplyr::mutate(newEnd=max(END)) %>%
#           dplyr::ungroup() %>%
#           dplyr::distinct(gid,.keep_all = T)  %>%
#           dplyr::mutate(START=newStart,END=newEnd) %>%
#           dplyr::select(-newEnd,-newStart)
#         
#         retain <- temp %>%
#           dplyr::filter(check==F & lag(check)==F)
#         
#         temp <- rbind(collapse,retain) %>%
#           dplyr::select(-gid,-check)
#         
#         j=j+1
#       }
#     }
#     print(head(temp))
#     all_collapsed[[i]] <- temp
#   }
#   return(all_collapsed)
# }
# 
# HDR_collapse_master <- getRegFreq(all_regions)
# 
# all_collapsed <- ldply(HDR_collapse_master, data.frame) %>%
#   dplyr::select(-strain) 
# 
# 
# # Adding resolution if genes are found in a HDR or not
# ortho_count_wCoord_HDR <- ortho_count_wCoord %>% 
#   dplyr::rowwise() %>%
#   dplyr::mutate(in_HDR = any(
#     contig == all_collapsed$CHROM &
#       start >= all_collapsed$START &
#       end <= all_collapsed$END
#   )) 
# 
# 
# stats <- ortho_count_wCoord_HDR %>%
#   dplyr::group_by(class) %>%
#   dplyr::summarise(count = n()) %>%
#   dplyr::ungroup()
# 
# ortho_count_wCoord_HDR_final <- ortho_count_wCoord_HDR %>%
#   dplyr::left_join(stats, by = 'class') %>%
#   dplyr::filter(contig != "MtDNA") %>%
#   dplyr::mutate(class = dplyr::recode(class, "core" = "core (14,799)", "accessory" = "accessory (4,660)", "private" = "private (462)")) %>%
#   dplyr::rename(Class = class, HDR = in_HDR)
# 
# 
# N2_genes_plot <- ggplot(ortho_count_wCoord_HDR_final) +
#   geom_point(aes(x = start / 1e6, y = freq * 100, color = Class, shape = HDR), size = 2) +
#   scale_color_manual(
#     name = "Gene set",
#     values = c("core (14,799)" = "green4", "accessory (4,660)" = "#DB6333", "private (462)" = "magenta3"),
#     limits = c("core (14,799)", "accessory (4,660)", "private (462)")) +
#   scale_shape_manual(
#     name = "In a HDR?",
#     values = c("TRUE" = 4, "FALSE" = 1)) +
#   geom_smooth(aes(x = start / 1e6, y = freq * 100), method = "loess", se = TRUE, color = "blue") +
#   facet_wrap(~contig, ncol = 6, scales = "free_x") +
#   labs(x = "Genome position (Mb)", y = "N2 gene frequency") +
#   scale_y_continuous(expand = c(0.01, 0), labels = scales::percent_format(scale = 1)) +
#   theme(
#     legend.position = 'right',
#     legend.text = element_text(size = 14, color = 'black'),
#     legend.title = element_text(size = 16, color = 'black'),
#     axis.text = element_text(size = 12, color = 'black'),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     panel.background = element_blank(),
#     panel.border = element_rect(fill = NA),
#     strip.text = element_text(size = 16, color = "black")
#   )
# N2_genes_plot
# 
# core_hdr <- ortho_count_wCoord_HDR_final %>%
#   dplyr::filter(Class == "core (14,799)") %>%
#   dplyr::count(HDR) # 1958/12841:InHDRs/notInHDRs - 13%
# 
# acc_hdr <- ortho_count_wCoord_HDR_final %>%
#   dplyr::filter(Class == "accessory (4,660)") %>%
#   dplyr::count(HDR) # 2273/2380:InHDRs/notInHDRs - 49%
# 
# priv_hdr <- ortho_count_wCoord_HDR_final %>%
#   dplyr::filter(Class == "private (462)") %>%
#   dplyr::count(HDR) # 105/357:InHDRs/notInHDRs - 23%
# 
# # ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/N2_geneFreq.png", N2_genes_plot, height = 10, width = 15, dpi = 600)
# 
# 
# ortho_count_wCoord_HDR_final$HDR <- factor(ortho_count_wCoord_HDR_final$HDR, levels = c(TRUE, FALSE))
# 
# ggplot(ortho_count_wCoord_HDR_final) +
#   geom_jitter(aes(x = HDR, y = freq * 100), width = 0.3, alpha = 0.5, size = 1) +
#   geom_boxplot(aes(x = HDR, y = freq * 100, fill = HDR), outlier.alpha = 0, alpha = 0.5) +
#   scale_fill_manual(
#     name = "Hyper-divergent region?",  
#     labels = c("non-HDR", "HDR"),  
#     values = c("#6495ED","#FF6347")  
#   ) +
#   labs(x="", y = "Frequency of all N2 genes in a HOG") +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) +
#   theme(
#         axis.text = element_text(size = 14, color = "black"),
#         axis.title = element_text(size = 16, color = "black"),
#         panel.grid = element_blank(),
#         legend.position = 'none',
#         legend.text = element_text(size = 8),
#         panel.background = element_blank(),
#         panel.border = element_rect(fill = NA))
#         # legend.key.height = unit(0.3, "cm"),  
#         # legend.key.width = unit(0.3, "cm"))





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
  dplyr::arrange(desc(count)) %>%
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
# mm_model <- nls(
#   num_priv_genes ~ Vmax * num_strains / (Km + num_strains),
#   data = rarefaction,
#   start = list(Vmax = max(rarefaction$num_priv_genes), Km = median(rarefaction$num_strains))
# )
# 
# # Extract parameters
# Vmax_est <- coef(mm_model)["Vmax"]
# Km_est <- coef(mm_model)["Km"]
# 
# # Define saturation proportion
# p <- 0.95  # 95% saturation
# 
# # Analytical saturation estimate
# x_sat <- (p / (1 - p)) * Km_est
# threshold <- p * Vmax_est
# additional_needed <- x_sat - max(rarefaction$num_strains)
# 
# # Predict over a reasonable x range for plotting
# x_pred <- seq(min(rarefaction$num_strains), x_sat * 1.1, length.out = 1000)
# y_pred <- predict(mm_model, newdata = data.frame(num_strains = x_pred))
# pred_df <- data.frame(num_strains = x_pred, num_priv_genes = y_pred)
# 
# # Print summary
# # cat("Estimated Vmax:", round(Vmax_est, 2), "\n")
# # cat("Estimated Km:", round(Km_est, 2), "\n")
# # cat("95% saturation reached at ~", round(x_sat, 2), "strains\n")
# # cat("Additional strains needed:", round(additional_needed, 2), "\n")
# 
# rfc_genes_fit <- ggplot(rarefaction, aes(x = num_strains, y = num_priv_genes)) +
#   geom_line(data = pred_df, aes(x = num_strains, y = num_priv_genes), color = "blue", size = 1) +
#   geom_point(color = "magenta3", size = 2, alpha = 0.8) +
#   geom_vline(xintercept = x_sat, linetype = "dashed", color = "black") +
#   geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
#   annotate("text", x = x_sat * 0.95, y = threshold,
#            label = paste0("Saturation at ~", round(x_sat), " strains"),
#            hjust = 0.9, vjust = 3, color = "blue", size = 4) +
#   xlab("Strains") +
#   ylab("Private genes") +
#   theme(
#     panel.background = element_blank(),
#     axis.title = element_text(size=13, color = 'black', face = 'bold'),
#     axis.text = element_text(size=11, color = 'black'),
#     panel.border = element_rect(fill = NA),
#     plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  
#   )
# rfc_genes_fit

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/rarefaction_115_labelled_fit.png", height = 5, width = 9, rfc_genes_fit, dpi = 600)









# ======================================================================================================================================================================================== #
# Plotting distribution of strains and count of private genes #
# ======================================================================================================================================================================================== #
trp <- private_final %>%
  dplyr::select(-Orthogroup) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("strain") %>%
  dplyr::mutate(strain = gsub("_count","", strain)) %>%
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

geneCount_OGs <- trp %>%
  dplyr::left_join(PC_count, by = "strain")

geneCount_OGs$strain <- factor(geneCount_OGs$strain, levels = geneCount_OGs$strain)
  
geneCount_OGs <- geneCount_OGs %>%
  dplyr::mutate(n_genes_scaled = n_genes / 100)

# mt_contigs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/nucmer_aln_WSs/numberOfContigsAln_MtDNA.tsv", col_names = F) %>%
#   tidyr::separate(X1, into = c("number","strain"), sep = ' ')
# 
# geneCount_OGs <- geneCount_OGs %>% dplyr::left_join(mt_contigs, by = 'strain') %>% dplyr::mutate(number = ifelse(is.na(number),0,number)) %>% dplyr::rename(mt_contigs = number) %>% 
#   dplyr::mutate(mt_contigs = ifelse(strain == "N2", 1, mt_contigs))
# 
# geneCount_OGs$strain <- factor(geneCount_OGs$strain, levels = geneCount_OGs$strain)

gene_mean <- geneCount_OGs %>% dplyr::summarise(mean_gene_count = mean(n_genes))

OGs_PCgenes <- ggplot(geneCount_OGs, aes(x = strain)) +
  geom_col(aes(y = count), fill = "magenta3", width = 0.6, alpha = 0.5, color = "black", linewidth = 0.3) +
  geom_line(aes(y = n_genes_scaled, group = 1), color = "blue", size = 1.2) +
  geom_point(aes(y = n_genes_scaled), color = "blue", size = 2) +
  geom_text(data = gene_mean, aes(x = 125, y = 250, label = paste0("Mean gene count: ", round(mean_gene_count))), size = 8, color = 'blue') +
  scale_y_continuous(
    expand = c(0.01, 0),
    name = "Private genes",
    sec.axis = sec_axis(~ ., name = "Protein-coding genes (1e02)")  
  ) +
  theme(
    axis.text.x = element_text(size = 8.5, color = 'black', angle = 70, hjust = 1),
    panel.background = element_blank(),
    axis.title = element_text(size = 22, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 18, color = 'black'),
    panel.border = element_rect(fill = NA),
    legend.position = "inside",
    legend.position.inside = c(0.9,0.8)
    # plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")  # 20pt on all sides
  ) +
  labs(x = "Strains")
OGs_PCgenes

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/HOGs_PCgenes.png", HOGs_PCgenes, height = 6, width = 10, dpi = 600)



corr <- geneCount_OGs %>% dplyr::select(strain, count) %>% dplyr::left_join(mt_contigs, by = 'strain') %>% dplyr::mutate(number = ifelse(is.na(number),0,number)) %>% dplyr::rename(mt_contigs = number, number_private = count) %>% dplyr::filter(strain != "N2") %>%
  dplyr::group_by(mt_contigs) %>%
  dplyr::mutate(mean_private = mean(number_private)) %>%
  dplyr::ungroup()

ggplot(corr) + 
  geom_point(aes(x = mt_contigs, y = number_private, color = strain), size = 2.5) + 
  geom_point(aes(x = mt_contigs, y = mean_private), color = 'black', size = 5) +
  geom_line(aes(x = mt_contigs, y = mean_private, group = 1), color = 'black') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black', face = 'bold')) +
  xlab("# of WS contigs that align to N2 MtDNA") +
  ylab("# of private genes")

# lm_model <- lm(count ~ n_genes_scaled, data = geneCount_HOGs)
# r2 <- summary(lm_model)$r.squared
# 
# blah <- ggplot(geneCount_HOGs, aes(x = n_genes, y = count)) +
#   geom_point(color = 'magenta3') +
#   geom_text(data = subset(geneCount_HOGs, n_genes > 27000 | n_genes < 20000 | count > 700), aes(label = strain), hjust = 0.9, vjust = 1.5, size = 5, color = "magenta3") +
#   geom_smooth(method = "lm", se = FALSE, color = "blue") + 
#   annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(r2, 3)), hjust = 2, vjust = 10, size = 7, color = "blue") +
#   labs(y = "Private genes", x = "Protein-coding genes") +
#   theme(
#     panel.background = element_blank(),
#     axis.title = element_text(size = 22, color = 'black', face = 'bold'),
#     axis.text = element_text(size = 18, color = 'black'),
#     panel.border = element_rect(fill = NA))
# blah

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/private_genes_PC.png", blah, height = 6, width = 10, dpi = 600)



# ======================================================================================================================================================================================== #
# Looking at relationship between number of private genes and number of single-exon genes #
# ======================================================================================================================================================================================== #

allGffs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/140_exons_strains.tsv", col_names = c("exonID", "strain", 'blah')) %>% dplyr::select(exonID,strain)

singleEx <- allGffs %>%
  dplyr::group_by(exonID, strain) %>%
  dplyr::mutate(num_exons = n()) %>%
  dplyr::filter(num_exons == 1) %>%
  dplyr::ungroup()

# awk '$3 == "exon" {print $NF}' c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3 | awk -F';' '{print $NF}' | sed 's|Parent=||'
N2_exons <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/N2_longestIso_exons_fixed.tsv", col_names = c("exonID", "strain")) %>%
  dplyr::group_by(exonID,strain) %>%
  dplyr::mutate(num_exons = n()) %>%
  dplyr::filter(num_exons == 1) %>%
  dplyr::ungroup() # 703 single-exon genes in N2

single_exonGenes <- singleEx %>%
  dplyr::bind_rows(N2_exons) %>%
  dplyr::count(strain, name = "num_SE_genes") %>%
  dplyr::filter(strain != "c_elegans")
  
priv_genes <- geneCount_OGs %>%
  dplyr::select(strain, count, n_genes) %>%
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
  # geom_text(data = subset(priv_genes, num_SE_genes > 3600 | num_SE_genes < 800), aes(x = privates, y = num_SE_genes, label = strain), hjust = 1.1, vjust = -0.5, size = 5, color = "firebrick") +
  # annotate("text", x = 1000, y = 1500, label = "PC genes", hjust = -1.5, vjust = -5, size = 5, color = "magenta") +
  # annotate("text", x = 1200, y = 4200, label = paste0("R² = ", round(r_sq, 3)), size = 5, color = "black") +
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
  dplyr::select(Orthogroup) %>%
  dplyr::pull()


# N2 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/N2_tran_gene.tsv", col_names = c('tran','gene'))

strain_SE <- singleEx %>%
  dplyr::bind_rows(N2_exons) %>%
  dplyr::select(-num_exons) %>%
  dplyr::mutate(exonID = gsub("\\.t[0-9]+", "", exonID)) %>%
  dplyr::mutate(exonID = gsub("transcript:", "", exonID)) %>%
  dplyr::mutate(exonID = ifelse(exonID %in% N2_tranGene$tran, N2_tranGene$gene, exonID)) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_single_exon = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(exonID = gsub("Parent=","",exonID)) %>%
  dplyr::select(strain, exonID, num_single_exon) 

privs <- ortho_genes_dd %>% 
  dplyr::bind_rows(private_OGs) %>%
  dplyr::filter(Orthogroup %in% private_ogs) %>%
  pivot_longer(
    cols = -Orthogroup,
    names_to = "strain",
    values_to = "genes"
  ) %>%
  dplyr::filter(!is.na(genes)) %>%                     
  dplyr::mutate(genes = str_split(genes, ",")) %>%     
  unnest(genes) %>%                             
  dplyr::arrange(strain) %>%
  dplyr::select(-Orthogroup) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_privs = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(genes = gsub("transcript_", "", genes)) %>% 
  dplyr::mutate(genes = stringr::str_trim(genes)) %>%
  dplyr::mutate(genes = ifelse(genes %in% N2_tranGene$tran, N2_tranGene$gene, genes)) %>%
  dplyr::mutate(genes = ifelse(strain != "N2", sub("\\..*","",genes), genes))
  
N2_privates <- privs %>% dplyr::filter(strain == "N2") %>% dplyr::select(genes)
# write.table(N2_privates,"/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/N2_private_genes.tsv", quote = F, col.names = F, row.names = F)

prop_SE_priv <- privs %>%
  inner_join(strain_SE, by = c("strain" = "strain", "genes" = "exonID")) %>%  
  group_by(strain) %>%
  summarise(
    num_priv_single_exon = n(),                        
    num_privs = first(num_privs),                      
    proportion = num_priv_single_exon / num_privs      
  )

okay2 <- prop_SE_priv %>% dplyr::filter(strain == "N2") # no private N2 genes are sinlge-exon genes

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
    axis.text.x = element_text(size = 16, color = 'black'),
    axis.title.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "right",
    panel.border = element_rect(fill = NA)) +
  scale_x_continuous(expand = c(0,0))
SE_plot























# ======================================================================================================================================================================================== #
# IPR analysis on genes that are private to N2 #
# ======================================================================================================================================================================================== #

# N2_tranGene <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/May_115_longestIso/N2_tran_gene.tsv", col_names = c('tran','gene')) %>%
  # dplyr::mutate(tran = paste0("transcript:", tran))

N2_tran_gene <- N2_tranGene %>%
  dplyr::mutate(tran = paste0("transcript:", tran))

N2_privates <- privs %>% dplyr::filter(strain == "N2") %>% dplyr::select(genes) %>% dplyr::pull()

ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/output/N2_IPR_allApps_20251019.tsv",  col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::left_join(N2_tran_gene, by = 'tran') %>%
  dplyr::select(-tran) %>%
  dplyr::select(gene,signature_accession,signature_description,IPR_accession,IPR_description,GO) %>%
  dplyr::rename(N2 = gene)

ipr_gene <- ipr %>%
  dplyr::filter(!is.na(IPR_description) & IPR_description != "-") %>%
  dplyr::select(N2, IPR_accession, IPR_description) %>%
  dplyr::distinct(N2, IPR_accession, IPR_description) # 14,996



##################################################################################################################
ipr_count <- ipr_gene %>% 
  dplyr::distinct(N2,IPR_accession, IPR_description) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(n_IPR_acc = n()) %>% 
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc)) %>%
  dplyr::select(-N2) %>% 
  dplyr::distinct() %>%
  dplyr::slice_head(n=30) 

ipr_count$IPR_description <- factor(ipr_count$IPR_description,levels = ipr_count$IPR_description[order(ipr_count$n_IPR_acc, decreasing = FALSE)])

ggplot(data = ipr_count, aes(x = n_IPR_acc, y = IPR_description)) +
  geom_bar(stat = "identity", fill = "magenta3", alpha = 0.5, color = "black", linewidth = 0.3) +
  theme(
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_blank(),
    legend.position = 'none',
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 18, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    plot.title = element_text(size = 20, color ='black', face = 'bold', hjust = 0.5),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)) +
  xlab("Number of genes with IPR term") +
  ggtitle("IPR terms for 531 private N2 genes") +
  scale_x_continuous(expand = c(0.01,0))
##################################################################################################################


# Define universe & HDR membership (annotated-only universe) 
univ_genes <- unique(ipr_gene$N2)
private_genes  <- intersect(N2_privates, univ_genes) # 251 / 531 genes.....

N <- length(univ_genes)
n <- length(private_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- ipr_gene %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- ipr_gene %>%
  dplyr::filter(N2 %in% private_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- ipr_gene %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig <- ipr_enrichment %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed <- ipr_gene %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(N2[N2 %in% private_genes]),
    genes_HDR   = paste(sort(unique(N2[N2 %in% private_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(N2),
    genes_all   = paste(sort(unique(N2)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "hyper-divergent regions")


# # Now for non-HDR arm genes
# univ_genes2 <- unique(ipr_gene$N2)
# nhdr_genes  <- intersect(nHD_gene_vector, univ_genes2) # 3,427 genes 
# 
# N <- length(univ_genes2)
# n <- length(nhdr_genes)
# 
# # Counts per IPR (k = in universe, x = in HDR subset) 
# k_tbl2 <- ipr_gene %>%
#   dplyr::count(IPR_accession, name = "k")
# 
# x_tbl2 <- ipr_gene %>%
#   dplyr::filter(N2 %in% nhdr_genes) %>%
#   dplyr::count(IPR_accession, name = "x")
# 
# desc_tbl2 <- ipr_gene %>%
#   dplyr::distinct(IPR_accession, IPR_description)
# 
# # Hypergeometric enrichment (one-sided)
# ipr_enrichment_nHDR <- k_tbl2 %>%
#   dplyr::left_join(x_tbl2, by = "IPR_accession") %>%
#   dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
#   dplyr::mutate(
#     pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
#     expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
#     enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
#     # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
#     OR = {
#       a <- x + 0.5                                # HDR & has IPR
#       b <- (n - x) + 0.5                          # HDR & no IPR
#       c <- (k - x) + 0.5                          # non-HDR & has IPR
#       d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
#       (a / b) / (c / d)
#     },
#     FDR_p.adjust = stats::p.adjust(pval, method = "BH")
#   ) %>%
#   dplyr::left_join(desc_tbl2, by = "IPR_accession") %>%
#   dplyr::mutate(N = N, n = n) %>%
#   dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
#   dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))
# 
# ipr_sig2 <- ipr_enrichment_nHDR %>%
#   dplyr::filter(FDR_p.adjust < 0.05)
# 
# ipr_sig2 %>% dplyr::slice_head(n = 20)
# 
# ipr_sig_gene_collapsed2 <- ipr_gene %>%
#   dplyr::filter(IPR_accession %in% ipr_sig2$IPR_accession) %>%
#   dplyr::group_by(IPR_accession) %>%
#   dplyr::summarise(
#     IPR_description = dplyr::first(stats::na.omit(IPR_description)),
#     n_genes_HDR = dplyr::n_distinct(N2[N2 %in% hdr_genes]),
#     genes_HDR   = paste(sort(unique(N2[N2 %in% hdr_genes])), collapse = ", "),
#     n_genes_all = dplyr::n_distinct(N2),
#     genes_all   = paste(sort(unique(N2)), collapse = ", "),
#     .groups = "drop") %>%
#   dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
#   dplyr::mutate(Region = "non HDRs")


binded <- ipr_sig_gene_collapsed %>% dplyr::arrange(FDR_p.adjust) 

data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=12, color='black'),
        axis.text.y = element_text(size=11, color='black'),
        # axis.title = element_text(size=12, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.76, 0.19),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 10, r = 10, l = 25, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  # shape = guide_legend(nrow=1, order = 3, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr






####################### GENE ONTOLOGY ##################################
go_ipr <- ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(N2, GO) %>% # 11,194 genes
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())

go_background <- go_ipr %>% dplyr::select(N2) %>% dplyr::distinct() %>% dplyr::pull()

# howmany <- go_ipr %>% dplyr::filter(N2 %in% N2_privates) %>% dplyr::distinct(N2) # 219 genes

GO_annotations <- AnnotationDbi::select(GO.db,
                                        keys=unique(go_ipr$GO),
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>%
  dplyr::rename(TERM = GOID, TERM_NAME = TERM)

merged_ont <- go_ipr %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) %>%
  dplyr::filter(!is.na(TERM_NAME))

# BP
enGO_HDR_merged_BP <- clusterProfiler::enricher(
  gene = N2_privates,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,N2),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = go_background,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_BP)

test <- as.data.table(enGO_HDR_merged_BP@result)

dotplot(enGO_HDR_merged_BP, showCategory = 40, title = "BP HDRs")

enGO_HDR_merged_plot_BP <- as.data.table(enGO_HDR_merged_BP@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(!is.na(Description)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(enGO_HDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=12, color='black'),
        axis.text.y = element_text(size=11, color='black'),
        axis.title = element_text(size=12, color='black', face = 'bold'),
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.76, 0.4),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = N2_privates,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,N2),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = go_background,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")

enGO_HDR_merged_plot <- as.data.table(enGO_HDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number()) %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (1)", Description))

plot_GO_MF <- ggplot(enGO_HDR_merged_plot) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.3),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF


final_plot <- cowplot::plot_grid(
  plot_ipr, plot_GO_BP,
  rel_heights = c(1,0.5),
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("a","b"),
  label_size = 14,
  label_fontface = "bold")
final_plot


