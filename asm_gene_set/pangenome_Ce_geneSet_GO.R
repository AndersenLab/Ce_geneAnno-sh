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
all_acc <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::filter(Orthogroup %in% acc_class) %>% tidyr::separate_rows(everything(), sep = ", ") %>% dplyr::select(-Orthogroup)
all_core <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::filter(Orthogroup %in% core_class) %>% tidyr::separate_rows(everything(), sep = ", ") %>% dplyr::select(-Orthogroup)

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
  dplyr::select(n_IPR_acc_background) %>%
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









#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - all pangenes as background

#==============================================================================================================================================================================================================================#
ipr_gene <- ipr %>%
  dplyr::filter(!is.na(IPR_description) & IPR_description != "-") %>%
  dplyr::select(QX1410, IPR_accession, IPR_description) %>%
  dplyr::distinct(QX1410,IPR_accession, IPR_description) %>% # 15,289
  dplyr::filter(QX1410 %in% arm_genes) # 5,301 genes!

# Define universe & HDR membership (annotated-only universe) 
univ_genes <- unique(ipr_gene$QX1410)
hdr_genes  <- intersect(HD_gene_vector, univ_genes) # 3,087 genes 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- ipr_gene %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- ipr_gene %>%
  dplyr::filter(QX1410 %in% hdr_genes) %>%
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
    n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
    genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(QX1410),
    genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "hyper-divergent regions")


# Now for non-HDR arm genes
univ_genes2 <- unique(ipr_gene$QX1410)
nhdr_genes  <- intersect(nHD_gene_vector, univ_genes2) # 3,427 genes 

N <- length(univ_genes2)
n <- length(nhdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl2 <- ipr_gene %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl2 <- ipr_gene %>%
  dplyr::filter(QX1410 %in% nhdr_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl2 <- ipr_gene %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment_nHDR <- k_tbl2 %>%
  dplyr::left_join(x_tbl2, by = "IPR_accession") %>%
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
  dplyr::left_join(desc_tbl2, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig2 <- ipr_enrichment_nHDR %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig2 %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed2 <- ipr_gene %>%
  dplyr::filter(IPR_accession %in% ipr_sig2$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
    genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(QX1410),
    genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "non HDRs")


binded <- ipr_sig_gene_collapsed %>% dplyr::bind_rows(ipr_sig_gene_collapsed2) %>% dplyr::arrange(FDR_p.adjust) 

data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size = 5, color='black', hjust = 1),
        legend.text = element_text(size = 4.5, color='black', hjust = 1),
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
  labs(title = "Enriched IPR terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr

#==============================================================================================================================================================================================================================#

# INTERPROSCAN (Gene Ontology) - all pangenes as background

#==============================================================================================================================================================================================================================#
go_ipr <- ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(QX1410, GO) %>% # 11,756 genes
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())


### Now with only arms as the background, not the entire genome
go_ipr_arms <- go_ipr %>% dplyr::filter(QX1410 %in% arm_genes)
IPR_GO_bckgrd_arms <- unique(go_ipr_arms$QX1410) # 3,905 genes


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
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
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
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=5, color='black', hjust = 1),
        legend.text = element_text(size=4.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.65, 0.42),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
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
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
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
        axis.title = element_text(size=10, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=5, color='black', hjust = 1),
        legend.text = element_text(size=4.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.65, 0.3),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(r = 10, b = 2, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF


final_plot <- cowplot::plot_grid(
  plot_ipr, plot_GO_BP, plot_GO_MF,
  rel_heights = c(1.5, 0.6, 1),
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("a","b","c"),
  label_size = 14,
  label_fontface = "bold")
final_plot
















# # ======================================================================================================================================================================================== #
# # Gene Ontology of N2 genes in each gene set #
# # ======================================================================================================================================================================================== #
# # Extract HOGs that have ONE N2 gene to lift over ontology to N2 orthologs #
# # all_relations_classification_rowid <- count %>% dplyr::rename(rowid = HOG)
# # all_relations_rowid <- all_relations %>% dplyr::rename(rowid = HOG)
# # ortho_genes_dd_rowid <- ortho_genes_dd %>% dplyr::rename(rowid = HOG)
# 
# genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
# all_genes_strain <- genes_strain %>%
#   dplyr::filter(strain != "N2" | grepl("protein_coding", attributes)) %>% #filter out non-protein coding N2 genes
#   dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
#   dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
#   dplyr::mutate(attributes = sub(";.*", "", attributes)) 
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
# write.table(genes_core, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/core.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 
# genes_acc <- n2_genes_GO_liftover %>%
#   dplyr::filter(class == "accessory") %>%
#   dplyr::pull(N2) # 4,236
# write.table(genes_acc, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/accessory.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 
# genes_private <- n2_genes_GO_liftover %>%
#   dplyr::filter(class == "private") %>%
#   dplyr::pull(N2) # 37
# write.table(genes_private, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/private.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# 
# 
# enrich_go <- function(wb_ids){
#   
#   mf <- enrichGO(gene          = wb_ids,
#                  OrgDb         = org.Ce.eg.db,
#                  ont           = "MF",
#                  pAdjustMethod = "bonferroni",
#                  keyType       = 'WORMBASE',
#                  pvalueCutoff  = 0.05,
#                  qvalueCutoff  = 0.05)
#   
#   bp <- enrichGO(gene          = wb_ids,
#                  OrgDb         = org.Ce.eg.db,
#                  ont           = "BP",
#                  pAdjustMethod = "bonferroni",
#                  keyType       = 'WORMBASE', 
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
#       scale_size_continuous(range = c(1,10), breaks = c(3.5)) + 
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
# 
# 
# 
# 
# 
# # ======================================================================================================================================================================================== #
# # RUNNING ENRICHMENT USING CLUSTERPROFILER ON EGGNOG-MAPPER OUTPUT  #
# # ======================================================================================================================================================================================== #
# tran_gene <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/longest_iso_tranName_WBGeneID.tsv", col_names = c("transcript","gene"))
# 
# ### Eukaryota ###
# N2_euk_bgd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output/euk_annotations/N2_longestIso_background_euk_GOannotations.tsv", col_names = c("transcript","GO", "term_name")) %>%
#   dplyr::left_join(tran_gene, by = "transcript") %>%
#   tidyr::separate_rows(GO, sep = ",") 
# 
# N2_euk_bgd_genes <- unique(N2_euk_bgd$gene) # a vector
# 
# GO_terms_euk <- AnnotationDbi::select(GO.db, keys=unique(N2_euk_bgd$GO), columns=c("TERM","DEFINITION"), keytype="GOID") %>% dplyr::rename(TERM = GOID, TERM_NAME = TERM)
# 
# # Load in gene set info #
# 
# core_test <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/core_tran_WBG.tsv", col_names = c("tran", "gene")) %>% 
#   dplyr::select(gene) %>% 
#   dplyr::pull() # need to vectorize
# 
# acc_test <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/accessory_tran_WBG.tsv", col_names = c("tran", "gene")) %>% 
#   dplyr::select(gene) %>%
#   dplyr::pull() # need to vectorize
# 
# priv_test <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/private_tran_WBG.tsv", col_names = c("tran", "gene")) %>% 
#   dplyr::select(gene) %>%
#   dplyr::pull() # need to vectorize
# 
# 
# # Core Eukaryota - No enrichment...?
# enGO_core_euk <- clusterProfiler::enricher(
#   gene = core_test,
#   TERM2GENE = N2_euk_bgd %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_euk %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_euk_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_core_euk)
# 
# dotplot(enGO_core_euk, showCategory = 40, title = "Core Eukaryota")
# 
# # Accessory Eukaryota
# enGO_acc_euk <- clusterProfiler::enricher(
#   gene = acc_test,
#   TERM2GENE = N2_euk_bgd %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_euk %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_euk_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_acc_euk)
# 
# dotplot(enGO_acc_euk, showCategory = 40, title = "Accessory Eukaryota")
# 
# # Private Eukaryota - none of the 37 private genes are found in the background set
# enGO_priv_euk <- clusterProfiler::enricher(
#   gene = priv_test,
#   TERM2GENE = N2_euk_bgd %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_euk %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_euk_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_priv_euk)
# 
# dotplot(enGO_priv_euk, showCategory = 40, title = "Private Eukaryota")
# 
# 
# 
# 
# 
# 
# ### Nematoda ###
# N2_nematoda_bgd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output/nematoda_annotations/N2_longestIso_background_nematoda_GOannotations.tsv", col_names = c("transcript","GO")) %>%
#   dplyr::left_join(tran_gene, by = "transcript") %>%
#   dplyr::select(gene,GO) %>%
#   tidyr::separate_rows(GO, sep = ",") 
# 
# N2_nematoda_bgd_genes <- unique(N2_nematoda_bgd$gene) # a vector
# 
# GO_terms_nematoda <- AnnotationDbi::select(GO.db, keys=unique(N2_nematoda_bgd$GO), columns=c("TERM","DEFINITION"), keytype="GOID") %>% dplyr::rename(TERM = GOID, TERM_NAME = TERM)
# 
# 
# # Core Nematoda 
# enGO_core_nematoda <- clusterProfiler::enricher(
#   gene = core_test,
#   TERM2GENE = N2_nematoda_bgd %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_nematoda %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_nematoda_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_core_nematoda)
# 
# dotplot(enGO_core_nematoda, showCategory = 40, title = "Core Nematoda")
# 
# 
# # Accessory Nematoda
# enGo_acc_nematoda <- clusterProfiler::enricher(
#   gene = acc_test,
#   TERM2GENE = N2_nematoda_bgd %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_nematoda %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_nematoda_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGo_acc_nematoda)
# 
# dotplot(enGo_acc_nematoda, showCategory = 40, title = "Accessory Nematoda")
# 
# 
# # Private Nematoda - none of the 37 genes in private have GO terms (not in the background set)
# enGo_priv_nematoda <- clusterProfiler::enricher(
#   gene = priv_test,
#   TERM2GENE = N2_nematoda_bgd %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_nematoda %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_nematoda_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGo_priv_nematoda)
# 
# dotplot(enGo_priv_nematoda, showCategory = 40, title = "Private Nematoda")
# 
# 
# # ======================================================================================================================================================================================== #
# # RUNNING ENRICHMENT USING CLUSTERPROFILER ON MERGED OUTPUT FROM EGGNOG-MAPPER AND ORG.CE.EG.DB #
# # ======================================================================================================================================================================================== #
# # I need to extract all WBGeneIDs and their GO annotation from org.Ce.ge.db
# # Which genes were identifed by eggNOG and not Ce.db, and vise versa, and which ones have no annotations
# gene_ids <- keys(org.Ce.eg.db, keytype = "WORMBASE")
# 
# n2_wbGenes <- tran_gene %>%
#   dplyr::select(gene) %>%
#   dplyr::distinct(gene)
# 
# # Documentation: https://www.bioconductor.org/packages/release/data/annotation/manuals/org.Ce.eg.db/man/org.Ce.eg.db.pdf
# gene2go <- AnnotationDbi::select(
#   org.Ce.eg.db,
#   keys = gene_ids,
#   columns = c("GO", "ONTOLOGY"),
#   keytype = "WORMBASE"
# )
# head(gene2go)
# 
# n2_gene_database <- gene2go %>%
#   dplyr::select(WORMBASE, GO, ONTOLOGY) %>%
#   dplyr::mutate(n2_gene = ifelse(WORMBASE %in% n2_wbGenes$gene, TRUE, FALSE)) %>%
#   dplyr::mutate(n2_gene = ifelse(n2_gene == TRUE, WORMBASE, n2_gene)) %>%
#   dplyr::mutate(eggNOG_GO = ifelse(n2_gene %in% N2_euk_bgd$gene, N2_euk_bgd$GO, NA))
# 
# 
# # n2_go_database <- n2_gene_database %>%
# #   dplyr::select(GO, n2_gene) %>%
# #   dplyr::filter(!is.na(GO)) # GO terms for 14,034 N2 genes
# 
# merged  <- n2_gene_database %>%
#   dplyr::select(n2_gene,GO,eggNOG_GO) %>%
#   dplyr::filter(n2_gene != FALSE) %>% 
#   dplyr::mutate(final_GO = dplyr::coalesce(GO, eggNOG_GO)) %>%
#   dplyr::select(n2_gene,final_GO) %>%
#   dplyr::filter(!is.na(final_GO)) %>%
#   dplyr::rename(gene = n2_gene, GO = final_GO) # 14,345 N2 genes
# 
# 
# N2_merged_bgd_genes <- unique(merged$gene) # a vector
# GO_terms_merged <- AnnotationDbi::select(GO.db, keys=unique(merged$GO), columns=c("TERM","DEFINITION"), keytype="GOID") %>% dplyr::rename(TERM = GOID, TERM_NAME = TERM)
# 
# ##^^^^^^ Can I extract ontology too? I.e., for only plotting BP, MF, or CC
# # Core merged
# enGO_core_merged <- clusterProfiler::enricher(
#   gene = core_test,
#   TERM2GENE = merged %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_merged %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_core_merged)
# 
# dotplot(enGO_core_merged, showCategory = 40, title = "Core Merged")
# 
# # Accessory merged
# enGO_acc_merged <- clusterProfiler::enricher(
#   gene = acc_test,
#   TERM2GENE = merged %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_merged %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_acc_merged)
# 
# dotplot(enGO_acc_merged, showCategory = 40, title = "Accessory Merged")
# 
# # Private merged
# enGO_priv_merged <- clusterProfiler::enricher(
#   gene = priv_test,
#   TERM2GENE = merged %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_merged %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_priv_merged)
# 
# dotplot(enGO_priv_merged, showCategory = 40, title = "Private Merged")
# 
# 
# 
