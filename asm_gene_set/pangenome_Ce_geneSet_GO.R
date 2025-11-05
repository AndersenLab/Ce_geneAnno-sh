library(readr)
library(org.Ce.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(clusterProfiler) ## BiocManager::install("clusterProfiler") # need this and the next package???
library(enrichplot)
library(cowplot)
library(GO.db)
library(AnnotationDbi)
library(stringr)
library(data.table)



# ======================================================================================================================================================================================== #
# OG matrix manipulation and plotting #
# ======================================================================================================================================================================================== #
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Oct16/Orthogroups/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))
# whatever <- ortho_genes_dd %>%
# dplyr::select(JU1581.braker.longest.protein, N2.longest.protein)

strainCol <- colnames(ortho_genes_dd)
ugh <- gsub(".20251012.inbred.blobFiltered.softMasked.braker.longestIso.protein","", strainCol)
ugh2 <- gsub(".20251014.inbred.blobFiltered.softMasked.braker.longestIso.protein","",ugh)
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


private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Oct16/Orthogroups/Orthogroups_UnassignedGenes.tsv") %>%
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




### Preparing lists of all core, accessory, and private pangenome genes
priv_class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::filter(class == "private") %>% dplyr::pull(Orthogroup)
################ Now just Hawaiian strains....
# HAWAII <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/Ce_speciesSheet_Hawaii.tsv", col_names = c("strain","desc","more_info")) %>% dplyr::mutate(strain = paste0(strain,"_count")) %>% dplyr::pull(strain)
# HAWAII_noCount <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/Ce_speciesSheet_Hawaii.tsv", col_names = c("strain","desc","more_info")) %>% dplyr::pull(strain)
# 
# private_freq <- 1/70
# priv_class <- all_relations %>%
#   dplyr::select(Orthogroup, 1, dplyr::any_of(HAWAII)) %>%
#   dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
#   dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
#   dplyr::mutate(freq = (sum / 70)) %>%
#   dplyr::mutate(
#     class = case_when(
#       freq == 1 ~ "core",
#       freq > private_freq & freq < 1 ~ "accessory",
#       freq == private_freq ~ "private",
#       TRUE ~ "undefined")) %>%
#   dplyr::filter(class == "private") %>% dplyr::pull(Orthogroup)
###########################################################
acc_class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::filter(class == "accessory") %>% dplyr::pull(Orthogroup)

core_class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::filter(class == "core") %>% dplyr::pull(Orthogroup)




all_priv <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::filter(Orthogroup %in% priv_class) %>% tidyr::separate_rows(everything(), sep = ", ") %>% dplyr::select(-Orthogroup)
# all_priv <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::filter(Orthogroup %in% priv_class) %>% dplyr::select(Orthogroup, 1, dplyr::any_of(HAWAII_noCount)) %>% tidyr::separate_rows(everything(), sep = ", ") %>% dplyr::select(-Orthogroup) # for HAWAII-only strains

### Private
## Adding "strain_" before each transcript to match IPR annotations
for (strain in names(all_priv)) {
  all_priv[[strain]] <- ifelse(is.na(all_priv[[strain]]), NA ,paste0(strain, "_", all_priv[[strain]]))
}

## Pulling all genes into a vector
all_priv_list <- c()
for (i in names(all_priv)) {
  vec <- all_priv[[i]]          # column as vector
  vec <- vec[!is.na(vec)]        # drop NAs
  all_priv_list <- c(all_priv_list, vec)  # append
}


### Accessory
## Adding "strain_" before each transcript to match IPR annotations
all_acc <- ortho_genes_dd %>% dplyr::filter(Orthogroup %in% acc_class) %>% 
  tidyr::pivot_longer(
    cols = -Orthogroup,
    names_to = "strain",
    values_to = "genes") %>%
  dplyr::filter(!is.na(genes)) %>%
  tidyr::separate_rows(genes, sep = ", ") %>%
  dplyr::mutate(genes = paste0(strain, "_", genes))

all_acc_list <- all_acc %>% dplyr::pull(genes)

### Core
## Adding "strain_" before each transcript to match IPR annotations
all_core <- ortho_genes_dd %>% dplyr::filter(Orthogroup %in% core_class) %>% 
  tidyr::pivot_longer(
      cols = -Orthogroup,
      names_to = "strain",
      values_to = "genes") %>%
      dplyr::filter(!is.na(genes)) %>%
      tidyr::separate_rows(genes, sep = ", ") %>%
      dplyr::mutate(genes = paste0(strain, "_", genes))
    
all_core_list <- all_core %>% dplyr::pull(genes)







# ======================================================================================================================================================================================== #
# IPR results
# ======================================================================================================================================================================================== #
onefourty_ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/pangenome/proteomes/output/140_iprAndGO.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways"))

n2_ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/output/N2_IPR_allApps_20251019.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::mutate(tran = gsub("transcript:","N2_transcript_",tran)) %>%
  dplyr::filter(IPR_accession != "-" | GO != "-")

all_ipr <- n2_ipr %>%
  dplyr::bind_rows(onefourty_ipr) %>%
  dplyr::select(tran, IPR_accession, IPR_description, GO)

all_ipr_background <- all_ipr %>% 
  dplyr::filter(!is.na(IPR_description) & IPR_accession != "-") %>%
  dplyr::distinct(tran, IPR_accession, IPR_description) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(n_IPR_acc_background = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_background))

ipr_background_genes <- all_ipr_background %>% dplyr::distinct(tran) %>% dplyr::pull()


# PRIVATE PANGENOME
priv_ipr_genes <- all_ipr_background %>%
  dplyr::filter(tran %in% all_priv_list) %>%
  dplyr::group_by(IPR_accession) %>% 
  dplyr::mutate(n_IPR_acc_priv = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_priv))

ratio <- all_ipr_background %>% dplyr::distinct(IPR_accession, .keep_all = T) %>%
  dplyr::left_join(priv_ipr_genes %>% dplyr::distinct(IPR_accession, .keep_all = T) %>% dplyr::select(IPR_accession, n_IPR_acc_priv), by = "IPR_accession") %>%
  dplyr::select(IPR_accession, IPR_description, n_IPR_acc_background, n_IPR_acc_priv) %>%
  dplyr::mutate(IPR_ratio = n_IPR_acc_priv / n_IPR_acc_background) %>%
  dplyr::filter(!is.na(n_IPR_acc_priv)) %>%
  dplyr::filter(n_IPR_acc_background > 100) %>% # removed very rarely appearing terms
  dplyr::arrange(desc(IPR_ratio))

ratio_plot <- ratio %>% dplyr::arrange(desc(IPR_ratio)) %>%
  dplyr::filter(IPR_accession != "IPR008164") %>% # filtering out "repeat of unknown function XGLTT
  dplyr::rename(`IPR count in private pangenome` = n_IPR_acc_priv) %>%
  dplyr::slice_head(n=50) %>%
  dplyr::mutate(IPR_description = factor(IPR_description, levels = rev(IPR_description)))

ggplot(data = ratio_plot, aes(x = IPR_ratio, y = IPR_description)) +
  geom_bar(aes(fill = `IPR count in private pangenome`),stat = "identity", alpha = 0.5, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = 'blue', high = 'magenta3') +
  theme(
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.7,0.2),
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 18, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    plot.title = element_text(size = 20, color ='black', face = 'bold', hjust = 0.5),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)) +
  xlab("Count of IPR in private / whole pangenome") +
  ggtitle("IPR enrichment for private pangenome") +
  scale_x_continuous(expand = c(0.01,0))

# ACCESSORY PANGENOME
acc_ipr_genes <- all_ipr_background %>%
  dplyr::filter(tran %in% all_acc_list) %>%
  dplyr::group_by(IPR_accession) %>% 
  dplyr::mutate(n_IPR_acc_acc = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_acc))

ratio <- all_ipr_background %>% dplyr::distinct(IPR_accession, .keep_all = T) %>%
  dplyr::left_join(acc_ipr_genes %>% dplyr::distinct(IPR_accession, .keep_all = T) %>% dplyr::select(IPR_accession, n_IPR_acc_acc), by = "IPR_accession") %>%
  dplyr::select(IPR_accession, IPR_description, n_IPR_acc_background, n_IPR_acc_acc) %>%
  dplyr::mutate(IPR_ratio = n_IPR_acc_acc / n_IPR_acc_background) %>%
  dplyr::filter(!is.na(n_IPR_acc_acc)) %>%
  # dplyr::filter(n_IPR_acc_background > 400) %>% # removed very rarely appearing terms
  dplyr::arrange(desc(IPR_ratio))

ratio_plot <- ratio %>% dplyr::arrange(desc(IPR_ratio)) %>%
  dplyr::filter(IPR_accession != "IPR008164") %>% # filtering out "repeat of unknown function XGLTT
  dplyr::rename(`IPR count in accessory pangenome` = n_IPR_acc_acc) %>%
  dplyr::slice_head(n=50) %>%
  dplyr::mutate(IPR_description = factor(IPR_description, levels = rev(IPR_description)))

ggplot(data = ratio_plot, aes(x = IPR_ratio, y = IPR_description)) +
  geom_bar(aes(fill = `IPR count in accessory pangenome`),stat = "identity", alpha = 0.5, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = 'blue', high = '#DB6333') +
  theme(
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.7,0.2),
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 18, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    plot.title = element_text(size = 20, color ='black', face = 'bold', hjust = 0.5),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)) +
  xlab("Count of IPR in accessory / whole pangenome") +
  ggtitle("IPR enrichment for accessory pangenome") +
  scale_x_continuous(expand = c(0.01,0))



# CORE PANGENOME
core_ipr_genes <- all_ipr_background %>%
  dplyr::filter(tran %in% all_core_list) %>%
  dplyr::group_by(IPR_accession) %>% 
  dplyr::mutate(n_IPR_acc_core = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_core))

ratio <- all_ipr_background %>% dplyr::distinct(IPR_accession, .keep_all = T) %>%
  dplyr::left_join(core_ipr_genes %>% dplyr::distinct(IPR_accession, .keep_all = T) %>% dplyr::select(IPR_accession, n_IPR_acc_core), by = "IPR_accession") %>%
  dplyr::select(IPR_accession, IPR_description, n_IPR_acc_background, n_IPR_acc_core) %>%
  dplyr::mutate(IPR_ratio = n_IPR_acc_core / n_IPR_acc_background) %>%
  dplyr::filter(!is.na(n_IPR_acc_core)) %>%
  dplyr::filter(n_IPR_acc_background > 100) %>% # removed very rarely appearing terms
  dplyr::arrange(desc(IPR_ratio))

ratio_plot <- ratio %>% dplyr::arrange(desc(IPR_ratio)) %>%
  dplyr::filter(IPR_accession != "IPR008164") %>% # filtering out "repeat of unknown function XGLTT
  dplyr::rename(`IPR count in core pangenome` = n_IPR_acc_core) %>%
  dplyr::slice_head(n=50) %>%
  dplyr::mutate(IPR_description = factor(IPR_description, levels = rev(IPR_description)))

ggplot(data = ratio_plot, aes(x = IPR_ratio, y = IPR_description)) +
  geom_bar(aes(fill = `IPR count in core pangenome`),stat = "identity", alpha = 0.5, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = 'blue', high = 'green4') +
  theme(
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.7,0.2),
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 18, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    plot.title = element_text(size = 20, color ='black', face = 'bold', hjust = 0.5),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)) +
  xlab("Count of IPR in core / whole pangenome") +
  ggtitle("IPR enrichment for core pangenome") +
  scale_x_continuous(expand = c(0.01,0))








#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - all pangenes as background

#==============================================================================================================================================================================================================================#
################################################# PRIVATE PANGENOME ##############################################################

# Define universe & HDR membership (annotated-only universe) 
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(all_priv_list, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% all_priv_list) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
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

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% all_priv_list]),
    genes_HDR   = paste(sort(unique(tran[tran %in% all_priv_list])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")

# 
# # Now for non-HDR arm genes
# univ_genes2 <- unique(ipr_gene$QX1410)
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
#   dplyr::filter(QX1410 %in% nhdr_genes) %>%
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
#     n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
#     genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
#     n_genes_all = dplyr::n_distinct(QX1410),
#     genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
#     .groups = "drop") %>%
#   dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
#   dplyr::mutate(Region = "non HDRs")


binded <- ipr_sig_gene_collapsed %>% dplyr::arrange(FDR_p.adjust) ##  %>% dplyr::bind_rows(ipr_sig_gene_collapsed2)

data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size = 7, color='black', hjust = 1),
        legend.text = element_text(size = 6, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 10, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  # shape = guide_legend(nrow=1, order = 3, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for private pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr





############################################################################### ACCESSORY PANGENOME #############################################################
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(all_acc_list, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k") %>%
  dplyr::mutate(k = as.numeric(k))

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% all_acc_list) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
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

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% all_acc_list]),
    genes_HDR   = paste(sort(unique(tran[tran %in% all_acc_list])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")

# 
# # Now for non-HDR arm genes
# univ_genes2 <- unique(ipr_gene$QX1410)
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
#   dplyr::filter(QX1410 %in% nhdr_genes) %>%
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
#     n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
#     genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
#     n_genes_all = dplyr::n_distinct(QX1410),
#     genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
#     .groups = "drop") %>%
#   dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
#   dplyr::mutate(Region = "non HDRs")


binded <- ipr_sig_gene_collapsed %>% dplyr::arrange(FDR_p.adjust) ##  %>% dplyr::bind_rows(ipr_sig_gene_collapsed2)

data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "#DB6333", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        # axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        # plot.title = element_blank(),
        legend.title = element_text(size = 7, color='black', hjust = 1),
        legend.text = element_text(size = 6, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 10, r = 10, l = 22, unit = "pt")) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  # shape = guide_legend(nrow=1, order = 3, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for accessory pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr








############################################################################### CORE PANGENOME #############################################################
univ_genes <- ipr_background_genes
hdr_genes  <- intersect(all_core_list, ipr_background_genes) 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- all_ipr_background %>%
  dplyr::count(IPR_accession, name = "k") %>%
  dplyr::mutate(k = as.numeric(k))

x_tbl <- all_ipr_background %>%
  dplyr::filter(tran %in% all_priv_list) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background %>%
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

ipr_sig_gene_collapsed <- all_ipr_background %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% all_core_list]),
    genes_HDR   = paste(sort(unique(tran[tran %in% all_core_list])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")

# 
# # Now for non-HDR arm genes
# univ_genes2 <- unique(ipr_gene$QX1410)
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
#   dplyr::filter(QX1410 %in% nhdr_genes) %>%
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
#     n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
#     genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
#     n_genes_all = dplyr::n_distinct(QX1410),
#     genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
#     .groups = "drop") %>%
#   dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
#   dplyr::mutate(Region = "non HDRs")


binded <- ipr_sig_gene_collapsed %>% dplyr::arrange(FDR_p.adjust) ##  %>% dplyr::bind_rows(ipr_sig_gene_collapsed2)

data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "green4", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size = 7, color='black', hjust = 1),
        legend.text = element_text(size = 6, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 10, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  # shape = guide_legend(nrow=1, order = 3, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for core pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr








#==============================================================================================================================================================================================================================#

# INTERPROSCAN (Gene Ontology) - all pangenes as background

#==============================================================================================================================================================================================================================#
go_ipr <- all_ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(tran, GO) %>% 
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())

### Now with only arms as the background, not the entire genome
IPR_GO_bckgrd_arms <- unique(go_ipr$tran) 

# how_many_HDR_GO_arm_genes <- go_ipr_arms %>% dplyr::filter(QX1410 %in% HD_gene_vector) # 2,105

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
  gene = all_priv_list,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,tran),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05)

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
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.x = element_text(size=14, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        # plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",s
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, t = 15, unit = "pt")) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = all_priv_list,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,tran),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
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
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.x = element_text(size=14, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        # plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, t = 15, unit = "pt")) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF




########################################### ACCESSORY PANGENOME GO ENRICHMENT ##########################################################
go_ipr <- all_ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(tran, GO) %>% 
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())

### Now with only arms as the background, not the entire genome
IPR_GO_bckgrd_arms <- unique(go_ipr$tran) 

# how_many_HDR_GO_arm_genes <- go_ipr_arms %>% dplyr::filter(QX1410 %in% HD_gene_vector) # 2,105

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
  gene = all_acc_list,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,tran),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05)

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
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(enGO_HDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "gold", high = "#DB6333", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.x = element_text(size=14, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        # plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.15),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, t = 15, unit = "pt")) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for accessory pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = all_acc_list,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,tran),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
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
  dplyr::slice_head(n = 30) #%>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (1)", Description))

plot_GO_MF <- ggplot(enGO_HDR_merged_plot) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.x = element_text(size=14, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        # plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, t = 15, unit = "pt")) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for accessory pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF





########################################### CORE PANGENOME GO ENRICHMENT ##########################################################
go_ipr <- all_ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(tran, GO) %>% 
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())

### Now with only arms as the background, not the entire genome
IPR_GO_bckgrd_arms <- unique(go_ipr$tran) 

# how_many_HDR_GO_arm_genes <- go_ipr_arms %>% dplyr::filter(QX1410 %in% HD_gene_vector) # 2,105

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
  gene = all_core_list,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,tran),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05)

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
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(enGO_HDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "yellowgreen", high = "green4", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.x = element_text(size=14, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        # plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, t = 15, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for core pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = all_core_list,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,tran),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
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
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.x = element_text(size=14, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        # plot.title = element_blank(),
        legend.title = element_text(size=14, color='black', hjust = 1),
        legend.text = element_text(size=12, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, t = 15, unit = "pt")) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for core pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF

