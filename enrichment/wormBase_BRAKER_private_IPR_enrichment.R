library(readr)
library(org.Ce.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(clusterProfiler) ## BiocManager::install("clusterProfiler") # need this and the next package???
library(enrichplot)
library(cowplot)

# ======================================================================================================================================================================================== #
# IPR results
# ======================================================================================================================================================================================== #
### N2 Wormbase
wormbase_ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/output/N2_IPR_allApps_20251019.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::select(tran, IPR_accession, IPR_description, GO) %>%
  dplyr::mutate(tran = gsub("transcript:","", tran))

wormbase_priv_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/withN2_BRAKER_orthofinder_run/N2_WormBase_specific_noN2BRAKER.tsv", col_names = 'gene') %>% 
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>%
  dplyr::pull()

all_ipr_background_wb <- wormbase_ipr %>% 
  dplyr::filter(!is.na(IPR_description) & IPR_accession != "-") %>%
  dplyr::distinct(tran, IPR_accession, IPR_description) %>% # only one type of IPR annotation per gene
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(n_IPR_acc_background = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_background)) %>%
  dplyr::mutate(tran = sub("\\.[^.]*$", "", tran)) # 14,996 genes have an IPR annotation

ipr_background_genes_wb <- all_ipr_background_wb %>% dplyr::distinct(tran) %>% dplyr::pull()

### N2 BRAKER
braker_ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/N2_BRAKER/output/N2_BRAKER_IPR_allApps_20260522.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::select(tran, IPR_accession, IPR_description, GO)

braker_priv_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/withN2_BRAKER_orthofinder_run/N2_BRAKER_specific_noN2wormBase.tsv", col_names = 'gene') %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>%
  dplyr::pull()

all_ipr_background_braker <- braker_ipr %>%
  dplyr::filter(!is.na(IPR_description) & IPR_accession != "-") %>%
  dplyr::distinct(tran, IPR_accession, IPR_description) %>% # only one type of IPR annotation per gene
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(n_IPR_acc_background = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_IPR_acc_background)) %>%
  dplyr::mutate(tran = sub("\\.[^.]*$", "", tran)) # 2,360,448 wild strain pangenome genes have an IPR annotation

ipr_background_genes_braker <- all_ipr_background_braker %>% dplyr::distinct(tran) %>% dplyr::pull()

#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - For WormBase genes in N2 that are missed by BRAKER

#==============================================================================================================================================================================================================================#
# Define universe & HDR membership (annotated-only universe)
univ_genes <- ipr_background_genes_wb
wb_priv_genes  <- intersect(wormbase_priv_genes, ipr_background_genes_wb)

N <- length(univ_genes) # 14996 / 19,984 genes received IPR annotations
n <- length(wb_priv_genes) # only 277 / 1,031 genes received IPR annotations - UNDERREPRESENTATION OF ANNOTATED GENES!

# Counts per IPR (k = in universe, x = in HDR subset)
k_tbl <- all_ipr_background_wb %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- all_ipr_background_wb %>%
  dplyr::filter(tran %in% wormbase_priv_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background_wb %>%
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

ipr_sig_gene_collapsed <- all_ipr_background_wb %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% wormbase_priv_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% wormbase_priv_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")


binded <- ipr_sig_gene_collapsed %>% dplyr::arrange(FDR_p.adjust) ##  %>% dplyr::bind_rows(ipr_sig_gene_collapsed2)

data_plt_wb <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr_wb <- ggplot(data_plt_wb) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt_wb$plotpoint, labels = data_plt_wb$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt_wb$n_genes_HDR, na.rm = TRUE)), round((max(data_plt_wb$n_genes_HDR, na.rm = TRUE) + min(data_plt_wb$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt_wb$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(4, 14), name = "Fold enrichment", breaks = pretty(data_plt_wb$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=14, color='black'),
        axis.text.y = element_text(size=16, color='black'),
        axis.title.x = element_text(size=18, color='black', face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_text(size = 18, color='black', hjust = 1),
        legend.text = element_text(size = 17, color='black', hjust = 1),
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
        panel.border = element_rect(fill = NA)) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  labs(title = "Enriched IPR terms for private pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr_wb

# The terms that appeared most in the N2 WormBase-specific genes:
counts_wb <- all_ipr_background_wb %>% dplyr::filter(tran %in% wormbase_priv_genes) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(count_in_priv = n()) %>%
  dplyr::ungroup() %>% dplyr::select(IPR_description, n_IPR_acc_background, count_in_priv) %>% dplyr::distinct()








#==============================================================================================================================================================================================================================#

# INTERPROSCAN enrichment test - For BRAKER genes that are predicted in N2 but not present in N2 WormBase annotations

#==============================================================================================================================================================================================================================#
# Define universe & HDR membership (annotated-only universe)
univ_genes <- ipr_background_genes_braker
wb_priv_genes  <- intersect(braker_priv_genes, ipr_background_genes_braker)

N <- length(univ_genes) # 16,748 / 22,479 genes received IPR annotations
n <- length(wb_priv_genes) # only 2,004 / 3,492 genes received IPR annotations 

# Counts per IPR (k = in universe, x = in HDR subset)
k_tbl <- all_ipr_background_braker %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- all_ipr_background_braker %>%
  dplyr::filter(tran %in% braker_priv_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- all_ipr_background_braker %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment_braker <- k_tbl %>%
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

ipr_sig_braker <- ipr_enrichment_braker %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig_braker %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed_braker <- all_ipr_background_braker %>%
  dplyr::filter(IPR_accession %in% ipr_sig_braker$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(tran[tran %in% braker_priv_genes]),
    genes_HDR   = paste(sort(unique(tran[tran %in% braker_priv_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(tran),
    genes_all   = paste(sort(unique(tran)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig_braker %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "Private pangenome")


binded_braker <- ipr_sig_gene_collapsed_braker %>% dplyr::arrange(FDR_p.adjust) ##  %>% dplyr::bind_rows(ipr_sig_gene_collapsed2)

data_plt_braker <- binded_braker %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr_braker <- ggplot(data_plt_braker) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt_braker$plotpoint, labels = data_plt_braker$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt_braker$n_genes_HDR, na.rm = TRUE)), round((max(data_plt_braker$n_genes_HDR, na.rm = TRUE) + min(data_plt_braker$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt_braker$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(4, 14), name = "Fold enrichment", breaks = pretty(data_plt_braker$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=14, color='black'),
        axis.text.y = element_text(size=16, color='black'),
        axis.title.x = element_text(size=18, color='black', face = 'bold'),
        plot.title = element_blank(),
        legend.title = element_text(size = 18, color='black', hjust = 1),
        legend.text = element_text(size = 17, color='black', hjust = 1),
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
        panel.border = element_rect(fill = NA)) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  labs(title = "Enriched IPR terms for private pangenome",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr_braker

# The terms that appeared most in the N2 WormBase-specific genes:
counts_braker <- all_ipr_background_braker %>% dplyr::filter(tran %in% braker_priv_genes) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::mutate(count_in_priv = n()) %>%
  dplyr::ungroup() %>% dplyr::select(IPR_description, n_IPR_acc_background, count_in_priv) %>% dplyr::distinct()




