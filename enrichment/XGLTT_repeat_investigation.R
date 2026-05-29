library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# ======================================================================================================================================================================================== #
# Load in IPR results
# ======================================================================================================================================================================================== #
all_ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/pangenome/proteomes/output/140WSs_andCGC1.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::select(tran, IPR_accession, IPR_description, GO)


xgltt_repeat <- all_ipr %>% dplyr::filter(IPR_description == "Repeat of unknown function XGLTT") %>% dplyr::distinct(tran,IPR_description) %>%
  tidyr::separate(tran, into = c("strain", "tran"), sep = "_") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_proteins = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain, num_proteins) %>%
  dplyr::arrange(desc(num_proteins)) 

order <- xgltt_repeat %>% dplyr::pull(strain)

geo_collect <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo.tsv") %>% dplyr::rename(strain = isotype)

geo.colors <- c("Hawaii"="#66C2A5", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="orange", "max_in_multiple_locations" = 'grey', 'Indian Ocean' = 'pink', "South America" = "cyan", "Asia" = "deeppink")

ggplot(xgltt_repeat %>% dplyr::left_join(geo_collect, by = "strain") %>% dplyr::mutate(geo = ifelse(strain == "CGC1", "Europe", geo)) %>% dplyr::mutate(strain = factor(strain, levels = order))) +
  geom_col(aes(x = strain, y = num_proteins, fill = geo), color = 'black') +
  scale_fill_manual(values = geo.colors) +
  theme(
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 13, color = 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 18, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.95,0.85),
    legend.text = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  labs(y = "Count", fill = "") +
  coord_cartesian(ylim = c(0,120))


# Gene set:
gene_set_res <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/ws_HDR_liftover/all_genes_class_OGs.tsv") %>% dplyr::select(strain,gene,orthogroup = Orthogroup,class)

xgltt_gene_set <- all_ipr %>% dplyr::filter(IPR_description == "Repeat of unknown function XGLTT") %>% dplyr::distinct(tran,IPR_description) %>%
  tidyr::separate(tran, into = c("strain", "gene"), sep = "_") %>% 
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>%
  dplyr::left_join(gene_set_res, by = c("strain","gene")) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(num_proteins = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(strain,class) %>%
  dplyr::mutate(prop_class = (n() / num_proteins) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,class,num_proteins,prop_class) %>%
  dplyr::mutate(class_prop_protein = num_proteins * (prop_class / 100)) %>%
  dplyr::mutate(strain = factor(strain, levels = order)) %>%
  dplyr::mutate(class = factor(class, levels = c("core", "private", "accessory")))

ggplot(xgltt_gene_set) +
  geom_col(aes(x = strain, y = class_prop_protein, fill = class), color = 'black', alpha = 0.5) +
  scale_fill_manual(values = c("core" = "green4", "accessory" = "#DB6333", "private" = "magenta3")) +
  theme(
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 13, color = 'black', angle = 60, hjust = 1),
    axis.title.y = element_text(size = 18, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'inside',
    legend.position.inside = c(0.95,0.85),
    legend.text = element_text(size = 16, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  labs(y = "Count", fill = "") +
  coord_cartesian(ylim = c(0,120))


core_ogs <- all_ipr %>% dplyr::filter(IPR_description == "Repeat of unknown function XGLTT") %>% dplyr::distinct(tran,IPR_description) %>%
  tidyr::separate(tran, into = c("strain", "gene"), sep = "_") %>% 
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene)) %>%
  dplyr::left_join(gene_set_res, by = c("strain","gene")) %>%
  dplyr::filter(class == 'core') %>% dplyr::distinct(orthogroup) %>% dplyr::pull()


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

core_n2 <- ortho_genes_dd %>% dplyr::filter(Orthogroup %in% core_ogs) %>%
  dplyr::select(N2) %>%
  tidyr::separate_rows(N2, sep = ", ") %>%
  dplyr::mutate(N2 = gsub("transcript_","", N2))

write.table(core_n2, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/pangenome/proteomes/output/XGLTT_N2_orthos_genes.tsv", col.names = F, row.names = F, quote = F)
  
  
  
