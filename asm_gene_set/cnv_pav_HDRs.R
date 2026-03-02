library(ggplot2)
library(dplyr)
library(data.table)
library(readr)

ws_hdr_ogs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/hdr_genes_OG_class.tsv")
strains <- ws_hdr_ogs %>% dplyr::distinct(strain) %>% dplyr::pull()
all_relations <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/OG_relations_matrix_count.tsv") %>%
  dplyr::rename_with(~ sub("_count$", "", .x))
  
  
cnv_pav_results_df = as.data.frame(matrix(ncol = 9, nrow = 140))
names(cnv_pav_results_df) = c("strain", "CNV_inHDR", "PAV_inHDR", "one_to_one_inHDR", "CNV_nonHDR", "PAV_nonHDR", "one_to_one_nonHDR", "HDR_OG_count", "nonHDR_OG_count")
  
for (i in 1:length(strains)) {
  soi <- strains[i]
  print(paste0("On strain: ", soi, ". ", i, "/140."))
  
  cnv_pav_results_df[i,1] = c(soi)
  
  hdr_ogs <- ws_hdr_ogs %>% dplyr::filter(strain == soi) %>% dplyr::distinct(Orthogroup) %>% dplyr::pull()
  # print(length(hdr_ogs))
  cnv_pav_results_df[i,8] = c(length(hdr_ogs))
  
  hdr_og_relations <- all_relations %>% dplyr::select(Orthogroup, N2, dplyr::all_of(soi)) %>% dplyr::filter(Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]]))) %>%
    dplyr::mutate(
      relation = dplyr::case_when(
        is.na(N2) | is.na(.data[[soi]]) ~ "PAV",
        N2 == .data[[soi]]              ~ "one_to_one",
        TRUE                            ~ "CNV" )) %>%
    dplyr::count(relation, name = "count")
  
  # print(hdr_og_relations %>% dplyr::distinct())
  
  cnv_hdr <- hdr_og_relations %>% dplyr::filter(relation == "CNV") %>% dplyr::pull(count)
  pav_hdr <- hdr_og_relations %>% dplyr::filter(relation == "PAV") %>% dplyr::pull(count)
  one_to_one <- hdr_og_relations %>% dplyr::filter(relation == "one_to_one") %>% dplyr::pull(count)
  
  cnv_pav_results_df[i,2] = c(cnv_hdr)
  cnv_pav_results_df[i,3] = c(pav_hdr)
  cnv_pav_results_df[i,4] = c(one_to_one)
  
  
  nonhdr_og_relations <- all_relations %>% dplyr::select(Orthogroup, N2, soi) %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]]))) %>%
    dplyr::mutate(
      relation = dplyr::case_when(
        is.na(N2) | is.na(.data[[soi]]) ~ "PAV",
        N2 == .data[[soi]]              ~ "one_to_one",
        TRUE                            ~ "CNV" )) %>%
    dplyr::count(relation, name = "count")
  
  nonhdr_ogs <- all_relations %>% dplyr::select(Orthogroup, N2, soi) %>% dplyr::filter(!Orthogroup %in% hdr_ogs) %>% filter(!(is.na(N2) & is.na(.data[[soi]])))
  
  cnv_pav_results_df[i,9] = c(nrow(nonhdr_ogs))
  
  cnv_hdr <- nonhdr_og_relations %>% dplyr::filter(relation == "CNV") %>% dplyr::pull(count)
  pav_hdr <- nonhdr_og_relations %>% dplyr::filter(relation == "PAV") %>% dplyr::pull(count)
  one_to_one <- nonhdr_og_relations %>% dplyr::filter(relation == "one_to_one") %>% dplyr::pull(count)
  
  cnv_pav_results_df[i,5] = c(cnv_hdr)
  cnv_pav_results_df[i,6] = c(pav_hdr)
  cnv_pav_results_df[i,7] = c(one_to_one)
}


# non-normalized results!

