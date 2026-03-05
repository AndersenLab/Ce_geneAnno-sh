library(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(rstatix)
library(ggpubr)
  
# Calculating stats for every wild strain and appending to a master dataframe
OG_enrichment <- function(ws_hdr_ogs, strains, og_matrix_relationships, single_copy_ogs) {
  
  cnv_pav_results_df = as.data.frame(matrix(ncol = 11, nrow = 140))
  names(cnv_pav_results_df) = c("strain", "CNV_inHDR", "PAV_inHDR", "one_to_one_inHDR", "single_copy_OGs_inHDRs", "CNV_nonHDR", "PAV_nonHDR", "one_to_one_nonHDR", "single_copy_OGs_nonHDRs", "HDR_OG_count", "nonHDR_OG_count")
  
  for (i in 1:length(strains)) {
    soi <- strains[i]
    print(paste0("On strain: ", soi, ". ", i, "/140."))
    
    cnv_pav_results_df[i,1] = c(soi)
    
    # HDRs stats
    hdr_ogs <- ws_hdr_ogs %>% dplyr::filter(strain == soi) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
    cnv_pav_results_df[i,10] = c(length(hdr_ogs))
  
    hdr_og_relations <- og_matrix_relationships %>% dplyr::select(Orthogroup, N2, dplyr::all_of(soi)) %>% dplyr::filter(Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]]))) %>%
      dplyr::mutate(
        relation = dplyr::case_when(
          is.na(N2) | is.na(.data[[soi]]) ~ "PAV",
          N2 == .data[[soi]]              ~ "one_to_one",
          TRUE                            ~ "CNV" )) %>%
      dplyr::count(relation, name = "count")
    
    sg_ogs <- og_matrix_relationships %>% dplyr::filter(Orthogroup %in% hdr_ogs) %>% dplyr::filter(Orthogroup %in% sc_ogs) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
  
    cnv_hdr <- hdr_og_relations %>% dplyr::filter(relation == "CNV") %>% dplyr::pull(count)
    pav_hdr <- hdr_og_relations %>% dplyr::filter(relation == "PAV") %>% dplyr::pull(count)
    one_to_one <- hdr_og_relations %>% dplyr::filter(relation == "one_to_one") %>% dplyr::pull(count)
    
    cnv_pav_results_df[i,2] = c(cnv_hdr)
    cnv_pav_results_df[i,3] = c(pav_hdr)
    cnv_pav_results_df[i,4] = c(one_to_one)
    cnv_pav_results_df[i,5] = c(length(sg_ogs))
    
    
    # outside of HDRs stats
    nonhdr_og_relations <- og_matrix_relationships %>% dplyr::select(Orthogroup, N2, soi) %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]]))) %>%
      dplyr::mutate(
        relation = dplyr::case_when(
          is.na(N2) | is.na(.data[[soi]]) ~ "PAV",
          N2 == .data[[soi]]              ~ "one_to_one",
          TRUE                            ~ "CNV" )) %>%
      dplyr::count(relation, name = "count")
    
    nonhdr_ogs <- og_matrix_relationships %>% dplyr::select(Orthogroup, N2, soi) %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]])))
    sg_ogs_nonHDR <- og_matrix_relationships %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% dplyr::filter(Orthogroup %in% single_copy_ogs) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
    
    
    cnv_pav_results_df[i,11] = c(nrow(nonhdr_ogs))
    
    cnv_hdr <- nonhdr_og_relations %>% dplyr::filter(relation == "CNV") %>% dplyr::pull(count)
    pav_hdr <- nonhdr_og_relations %>% dplyr::filter(relation == "PAV") %>% dplyr::pull(count)
    one_to_one <- nonhdr_og_relations %>% dplyr::filter(relation == "one_to_one") %>% dplyr::pull(count)
    
    cnv_pav_results_df[i,6] = c(cnv_hdr)
    cnv_pav_results_df[i,7] = c(pav_hdr)
    cnv_pav_results_df[i,8] = c(one_to_one)
    cnv_pav_results_df[i,9] = c(length(sg_ogs_nonHDR))
  }
  return(cnv_pav_results_df)
}

# Preparing data:
ws_hdr_ogs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/hdr_genes_OG_class.tsv")
strains <- ws_hdr_ogs %>% dplyr::distinct(strain) %>% dplyr::pull()
all_relations <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/OG_relations_matrix_count.tsv") %>%
  dplyr::rename_with(~ sub("_count$", "", .x))
sc_ogs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups_SingleCopyOrthologues.txt", col_names = "OGs") %>%
  dplyr::pull(OGs)

# Running orthogroup enrichment
og_enrich_results <- OG_enrichment(ws_hdr_ogs, strains, all_relations, sc_ogs)

### Diagnostic plots
# # non-normalized results plotted
# plot_df <- og_enrich_results %>%
#   tidyr::pivot_longer(
#     cols = -strain,
#     names_to = "metric",
#     values_to = "value") %>%
#   dplyr::mutate(region = case_when(
#       stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
#       stringr::str_detect(metric, "HDR") ~ "HDR",
#       TRUE ~ NA)) %>%
#   dplyr::filter(!is.na(region)) %>%
#   dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
#                 metric = ifelse(metric == "HDR_OG_count", "OG_count", 
#                                 ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
#                 stat = metric %>%
#                   stringr::str_remove("_inHDRs?$") %>%
#                   stringr::str_remove("_nonHDRs?$"),
#                 stat = ifelse(stat == "one_to_one","one-to-one", 
#                               ifelse(stat == "single_copy_OGs", "single-copy OGs",
#                                      ifelse(stat == "OG_count", "OGs", stat))),
#                 stat = factor(stat, levels = c("CNV","PAV","one-to-one","single-copy OGs","OGs")))
# 
# ggplot(plot_df, aes(x = stat, y = value, fill = region)) +
#   geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
#   scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 22, color = 'black', face = 'bold'),
#     panel.grid.major.x = element_blank(),
#     legend.title = element_blank(),
#     legend.text = element_text(size = 20, color = 'black'),
#     axis.text.y = element_text(size = 18, color = 'black'),
#     axis.title.y = element_text(size = 22, color = 'black', face = 'bold')
#   ) +
#   scale_y_log10() +
#   labs(y = "Orthogroup count", fill = "Region")
# 
# 
# # normalizing for # of orthogroups in HDRs versus not in HDRs
# plot_df_norm <- og_enrich_results %>%
#   dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
#   dplyr::mutate(
#     dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
#     dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
#   tidyr::pivot_longer(
#     cols = -strain,
#     names_to = "metric",
#     values_to = "value") %>%
#   dplyr::mutate(region = case_when(
#     stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
#     stringr::str_detect(metric, "HDR") ~ "HDR",
#     TRUE ~ NA)) %>%
#   dplyr::filter(!is.na(region)) %>%
#   dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
#                 metric = ifelse(metric == "inHDR_OG_count", "OG_count", 
#                                 ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
#                 stat = metric %>%
#                   stringr::str_remove("_inHDRs?$") %>%
#                   stringr::str_remove("_nonHDRs?$"),
#                 stat = ifelse(stat == "one_to_one","n-to-n", 
#                               ifelse(stat == "single_copy_OGs", "single-copy OGs",
#                                      ifelse(stat == "OG_count", "OGs", stat))),
#                 stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
#   dplyr::filter(stat != "OGs")
# 
# 
# ggplot(plot_df_norm, aes(x = stat, y = value, fill = region)) +
#   geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
#   scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
#   theme_bw() +
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 26, color = 'black', face = 'bold'),
#     panel.grid.major.x = element_blank(),
#     legend.box.background = element_rect(color = "black", size = 1),
#     legend.position  = 'inside',
#     legend.position.inside = c(0.945,0.945),
#     legend.title = element_blank(),
#     legend.key.size = unit(1.5, "cm"),
#     legend.text = element_text(size = 24, color = 'black'),
#     axis.text.y = element_text(size = 18, color = 'black'),
#     axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
#   ) +
#   scale_y_log10() +
#   labs(y = "Orthogroup count", fill = "Region")


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm <- og_enrich_results %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs")

wilcox_results <- plot_df_norm  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos <- plot_df_norm %>%
  group_by(stat) %>%
  summarise(y.position = max(value, na.rm = TRUE) * 1.3, .groups = "drop")

wilcox_results <- wilcox_results %>%
  left_join(y_pos, by = "stat")

ggplot(plot_df_norm, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Normalized orthogroup count", fill = "Region") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill =NA),
    legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  # scale_x_discrete(labels = c(
  #   "CNV" = "CNV",
  #   "PAV" = "PAV",
  #   "n-to-n" = expression(italic(n) * "-to-" * italic(n)),
  #   "single-copy OGs" = "single-copy OGs"
  # ))
  scale_x_discrete(labels = c(
  "CNV" = expression(bold("CNV")),
  "PAV" = expression(bold("PAV")),
  "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
  "single-copy OGs" = expression(bold("single-copy OGs"))
  ))





#####################################################################################################
# Adding on gene-class specific variation
#####################################################################################################
### GPCRs ###
ws_gpcrs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/140_IPR_gpcrs.tsv") %>% dplyr::mutate(gpcr = TRUE) # see line 1664 in orthogroup_vis.R

gpcr_ws_hdr_ogs <- ws_hdr_ogs %>% dplyr::left_join(ws_gpcrs, by = c("strain","gene")) %>% dplyr::filter(gpcr == TRUE)

all_ws_genes_class_og <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv")

all_ws_gpcr_OGs <- all_ws_genes_class_og %>% dplyr::left_join(ws_gpcrs, by = c("strain","gene")) %>% dplyr::filter(gpcr == TRUE) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()

gpcr_all_relations <- all_relations %>% dplyr::filter(Orthogroup %in% all_ws_gpcr_OGs)


# Running orthogroup enrichment
og_enrich_results_GPCRs <- OG_enrichment(gpcr_ws_hdr_ogs, strains, gpcr_all_relations, sc_ogs)


#### Calculating stats and adding to plot: 
# Paired Wilcoxon signed-rank test
plot_df_norm_gpcrs <- og_enrich_results_GPCRs %>%
  dplyr::rename(inHDR_OG_count = HDR_OG_count) %>%
  dplyr::mutate(
    dplyr::across(dplyr::matches("inHDR"), ~ .x / inHDR_OG_count),
    dplyr::across(dplyr::matches("nonHDR"), ~ .x / nonHDR_OG_count)) %>%
  tidyr::pivot_longer(
    cols = -strain,
    names_to = "metric",
    values_to = "value") %>%
  dplyr::mutate(region = case_when(
    stringr::str_detect(metric, "nonHDR") ~ "non-HDR",
    stringr::str_detect(metric, "HDR") ~ "HDR",
    TRUE ~ NA)) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::mutate(region = factor(region, levels = c("HDR", "non-HDR")),
                metric = ifelse(metric == "inHDR_OG_count", "OG_count",
                                ifelse(metric == "nonHDR_OG_count", "OG_count", metric)),
                stat = metric %>%
                  stringr::str_remove("_inHDRs?$") %>%
                  stringr::str_remove("_nonHDRs?$"),
                stat = ifelse(stat == "one_to_one","n-to-n",
                              ifelse(stat == "single_copy_OGs", "single-copy OGs",
                                     ifelse(stat == "OG_count", "OGs", stat))),
                stat = factor(stat, levels = c("CNV","PAV","n-to-n","single-copy OGs","OGs"))) %>%
  dplyr::filter(stat != "OGs")

wilcox_results <- plot_df_norm_gpcrs  %>%
  dplyr::group_by(stat) %>%
  wilcox_test(value ~ region, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>% # Benjamini-Hochberg correction
  add_significance()

y_pos <- plot_df_norm_gpcrs %>%
  group_by(stat) %>%
  summarise(y.position = max(value, na.rm = TRUE) * 1.3, .groups = "drop")

wilcox_results <- wilcox_results %>%
  left_join(y_pos, by = "stat")

ggplot(plot_df_norm_gpcrs, aes(x = stat, y = value, fill = region)) +
  geom_boxplot(outlier.size = 0.6, width = 0.7, position = position_dodge(width = 0.75), outlier.shape = NA, alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1.25) +
  scale_fill_manual(values = c("HDR" = "red", "non-HDR" = "blue")) +
  stat_pvalue_manual(wilcox_results, label = "p.adj.signif", x = "stat", y.position = "y.position", inherit.aes = FALSE, color = 'black', size = 8) +
  scale_y_log10() +
  labs(y = "Normalized orthogroup count", fill = "Region") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 26, color = 'black'),
    panel.grid.major.x = element_blank(),
    legend.box.background = element_rect(color = "black", size = 1),
    legend.position  = 'inside',
    panel.border = element_rect(color = 'black', fill =NA),
    legend.position.inside = c(0.945, 0.945),
    legend.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(size = 24, color = 'black'),
    axis.text.y = element_text(size = 18, color = 'black'),
    axis.title.y = element_text(size = 26, color = 'black', face = 'bold')
  )  +
  # scale_x_discrete(labels = c(
  #   "CNV" = "CNV",
  #   "PAV" = "PAV",
  #   "n-to-n" = expression(italic(n) * "-to-" * italic(n)),
  #   "single-copy OGs" = "single-copy OGs"
  # ))
  scale_x_discrete(labels = c(
    "CNV" = expression(bold("CNV")),
    "PAV" = expression(bold("PAV")),
    "n-to-n" = expression(bold(bolditalic(n) * "-to-" * bolditalic(n))),
    "single-copy OGs" = expression(bold("single-copy OGs"))
  ))



























