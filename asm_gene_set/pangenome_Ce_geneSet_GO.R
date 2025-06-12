library(readr)
library(org.Ce.eg.db)
library(dplyr)
library(ggplot2)
library(tidyr)
library(clusterProfiler) ## BiocManager::install("clusterProfiler") # need this and the next package???
library(enrichplot)
library(cowplot)
library(GO.db)


# ======================================================================================================================================================================================== #
# HOG matrix manipulation and plotting #
# ======================================================================================================================================================================================== #

# Converting transcripts to genes and removing all duplicate genes (not needed anymore, using longest isoform) in a dataframe cell:
# see script - /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/tran_gene.sh 
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0504_genes.tsv")

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

all_relations <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))


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
  dplyr::count(freq, class) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) 


gs_allOrtho <- ggplot(data = classification, aes(x = freq * 100, y = percent, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("HOGs") + # longest isoform
  xlab("Frequency") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    plot.margin = margin(l = 20),
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=16, color = 'black'),
    axis.text = element_text(size=14, color = 'black')
  )
gs_allOrtho

HOG_class_count <- classification %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_HOG = sum(n)) %>%
  dplyr::ungroup()




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
# ======================================================================================================================================================================================== #
# Gene Ontology of N2 genes in each gene set #
# ======================================================================================================================================================================================== #
# Extract HOGs that have ONE N2 gene to lift over ontology to N2 orthologs #
# all_relations_classification_rowid <- count %>% dplyr::rename(rowid = HOG)
# all_relations_rowid <- all_relations %>% dplyr::rename(rowid = HOG)
# ortho_genes_dd_rowid <- ortho_genes_dd %>% dplyr::rename(rowid = HOG)

# Extract rowids for HOGs that have one N2 gene contributing
oneN2 <- all_relations %>%
  dplyr::filter(N2_count == 1) %>% 
  dplyr::pull(HOG)

private_n2 <- count %>%
  dplyr::filter(class == "private") %>%
  dplyr::filter(!is.na(N2_count)) %>%
  dplyr::pull(HOG) # 16 - this is correct

one_n2_and_private <- c(oneN2,private_n2)

HOG_classification <- count %>%
  dplyr::select(HOG,freq,class)

n2_genes_GO_liftover <- ortho_genes_dd %>% dplyr::filter(HOG %in% one_n2_and_private) %>%
  dplyr::select(-OG, -"Gene Tree Parent Clade") %>%
  dplyr::left_join(HOG_classification, by = "HOG") %>%
  tidyr::separate_rows(N2, sep = ",\\s*") # need to split N2 genes in private HOGs because there are multiple N2 genes 

all_N2_genes <- all_genes_strain %>%
  dplyr::filter(strain == "N2") %>%
  dplyr::pull(attributes)

genes_core <- n2_genes_GO_liftover %>%
  dplyr::filter(class == "core") %>%
  dplyr::pull(N2) # 14,550
write.table(genes_core, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/core.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

genes_acc <- n2_genes_GO_liftover %>%
  dplyr::filter(class == "accessory") %>%
  dplyr::pull(N2) # 4,236
write.table(genes_acc, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/accessory.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

genes_private <- n2_genes_GO_liftover %>%
  dplyr::filter(class == "private") %>%
  dplyr::pull(N2) # 37
write.table(genes_private, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/private.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


enrich_go <- function(wb_ids){
  
  mf <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "bonferroni",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  bp <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "bonferroni",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  return(list(mf, bp))
}



N2_anno <- enrich_go(all_N2_genes)
core_anno <- enrich_go(genes_core)
acc_anno <- enrich_go(genes_acc)
private_anno <- enrich_go(genes_private)

data_man_plot <- function(N2_anno, geneSet_anno, gene_set) {
  genesBP_N2 <- setReadable(N2_anno[[2]], OrgDb = org.Ce.eg.db)
  gene_BP_geneSet <- setReadable(geneSet_anno[[2]], OrgDb = org.Ce.eg.db)
  
  df_GO_enrich_BP <- rbind(dplyr::mutate(genesBP_N2[], freq="control_N2_genes"), dplyr::mutate(gene_BP_geneSet[], freq = paste0(gene_set,"_genes"))) %>%
    tidyr::separate(GeneRatio, c("genes_enrich", "genes_in_database")) %>%
    tidyr::separate(BgRatio, c("genes_in_geneSet", "N2_genes_in_database")) %>%
    dplyr::mutate(genes_enrich = as.numeric(genes_enrich), genes_in_database = as.numeric(genes_in_database), genes_in_geneSet = as.numeric(genes_in_geneSet), N2_genes_in_database = as.numeric(N2_genes_in_database)) %>%
    dplyr::mutate(GeneRatio = genes_enrich/genes_in_database, BgRatio = genes_in_geneSet/N2_genes_in_database) %>%
    dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
    dplyr::arrange(p.adjust)
  
  GO_list_BP <- df_GO_enrich_BP$Description
  
  GO_list_BP_plotpoint <- data.frame(Description=GO_list_BP, plotpoint=length(GO_list_BP):1)
  
  df_GO_enrich_BP_sum <- df_GO_enrich_BP %>%
    dplyr::filter(Description %in% GO_list_BP) %>%
    dplyr::left_join(., GO_list_BP_plotpoint, by = "Description") %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(class_gene_total = max(genes_in_geneSet)) %>%
    dplyr::ungroup()
  
  df_class_total_BP <- df_GO_enrich_BP_sum %>%
    dplyr::distinct(Description, class_gene_total, plotpoint) %>%
    dplyr::arrange(-plotpoint)
  
  df_GO_enrich_BP_sum$freq <- factor(df_GO_enrich_BP_sum$freq, levels = c(paste0(gene_set,"_genes"),"control_N2_genes"), labels = c(gene_set,"All N2 genes"))
  
  if (gene_set == "Core") { 
    plot_GO_BP <- ggplot(df_GO_enrich_BP_sum) +
      geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
      geom_point(aes(x = -log10(p.adjust), y = plotpoint, shape = freq, size = EnrichRatio, fill = Count)) +
      scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.02,0.02)) +
      scale_fill_gradient(low = "darkolivegreen1", high = "green4", breaks = c(round(min(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(median(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(max(df_GO_enrich_BP_sum$Count, na.rm = TRUE)))) +
      scale_shape_manual(values=c(21, 22)) +
      scale_size_continuous(range = c(1,4)) + #breaks = c(1, 2, 3)) + 
      theme(axis.text.x = element_text(size=9, color='black'), 
            axis.text.y = element_text(size=7, color='black'),
            axis.title = element_text(size=10, color='black', face = 'bold'), 
            legend.title = element_text(size=7, color='black'), 
            legend.text = element_text(size=6, color='black'), 
            legend.position = "inside", 
            legend.position.inside = c(0.75, 0.18),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.direction = "horizontal", legend.box = "vertical",
            plot.margin = margin(l = 20, unit = "pt"),
            legend.spacing.y = unit(0.02, 'in'),
            # legend.spacing.x = unit(0.02, "cm"),
            legend.key.size = unit(0.3, "cm")) +
      guides(shape = guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3), force = TRUE), 
             size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE), 
             fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2, force = TRUE)) +
      labs(title = paste0(gene_set," gene set"), x="-log10 (corrected p-value)", shape = "Gene Set", size = "Fold enrichment", fill = "Gene counts") 
  } else if (gene_set == "Accessory") {
    plot_GO_BP <- ggplot(df_GO_enrich_BP_sum) +
      geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
      geom_point(aes(x = -log10(p.adjust), y = plotpoint, shape = freq, size = EnrichRatio, fill = Count)) +
      scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.02,0.02)) +
      scale_fill_gradient(low = "burlywood1", high = "#DB6333", breaks = c(round(min(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(median(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(max(df_GO_enrich_BP_sum$Count, na.rm = TRUE)))) +
      scale_shape_manual(values=c(21, 22)) +
      scale_size_continuous(range = c(1,10), breaks = c(3.5)) + 
      theme(axis.text.x = element_text(size=9, color='black'), 
            axis.text.y = element_text(size=7, color='black'),
            axis.title = element_text(size=10, color='black', face = 'bold'), 
            legend.title = element_text(size=7, color='black'), 
            legend.text = element_text(size=6, color='black'), 
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            
            legend.position = "inside", 
            legend.position.inside = c(0.72, 0.18),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.direction = "horizontal", legend.box = "vertical",
            # plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
            legend.spacing.y = unit(0.02, 'in'),
            # legend.spacing.x = unit(0.02, "cm"),
            legend.key.size = unit(0.3, 'cm')) +
      guides(shape = guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3), force = TRUE), 
             size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE), 
             fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2, force = TRUE)) +
      labs(title = paste0(gene_set," gene set"), x="-log10 (corrected p-value)", shape = "Gene Set", size = "Fold enrichment", fill = "Gene counts") 
  } else if (gene_set == "Private") {
    plot_GO_BP <- ggplot(df_GO_enrich_BP_sum) +
      geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
      geom_point(aes(x = -log10(p.adjust), y = plotpoint, shape = freq, size = EnrichRatio, fill = Count)) +
      scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.02,0.02)) +
      scale_fill_gradient(low = "thistle1", high = "magenta3", breaks = c(round(min(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(median(df_GO_enrich_BP_sum$Count, na.rm = TRUE)), round(max(df_GO_enrich_BP_sum$Count, na.rm = TRUE)))) +
      scale_shape_manual(values=c(21, 22)) +
      scale_size_continuous(range = c(2,8), breaks = c(50, 150, 250)) + 
      theme(axis.text.x = element_text(size=9, color='black'), 
            axis.text.y = element_text(size=7, color='black'),
            axis.title = element_text(size=10, color='black', face = 'bold'), 
            legend.title = element_text(size=7, color='black'), 
            legend.text = element_text(size=6, color='black'), 
            legend.position = "inside", 
            legend.position.inside = c(0.72, 0.18),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            legend.direction = "horizontal", legend.box = "vertical",
            # plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
            legend.spacing.y = unit(0.02, 'in'),
            #legend.spacing.x = unit(0.02, "cm"),
            legend.key.size = unit(0.3, 'cm')) +
      guides(shape = guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3), force = TRUE), 
             size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE), 
             fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2, force = TRUE)) +
      labs(title = paste0(gene_set," gene set"), x="-log10 (corrected p-value)", shape = "Gene Set", size = "Fold enrichment", fill = "Gene counts") 
  } else {
    stop("Please provide gene set as 'Core', 'Accessory', or 'Private'")
  }
  
  return(plot_GO_BP)
}

core_plot <- data_man_plot(N2_anno, core_anno, "Core")
core_plot

acc_plot <- data_man_plot(N2_anno, acc_anno, "Accessory")
acc_plot

private_plot <- data_man_plot(N2_anno, private_anno, "Private")
private_plot



GO_all <- plot_grid(
  core_plot, acc_plot, private_plot,
  nrow = 1,
  rel_heights = c(1, 1), rel_widths = c(1,1,1),
  align = "hv"
)
GO_all


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/GO_pangenomeGeneSet.png", GO_all, height = 8, width = 26, dpi = 600)







# ======================================================================================================================================================================================== #
# RUNNING ENRICHMENT USING CLUSTERPROFILER ON EGGNOG-MAPPER OUTPUT  #
# ======================================================================================================================================================================================== #
tran_gene <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/longest_iso_tranName_WBGeneID.tsv", col_names = c("transcript","gene"))

### Eukaryota ###
N2_euk_bgd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output/euk_annotations/N2_longestIso_background_euk_GOannotations.tsv", col_names = c("transcript","GO", "term_name")) %>%
  dplyr::left_join(tran_gene, by = "transcript") %>%
  tidyr::separate_rows(GO, sep = ",") 

N2_euk_bgd_genes <- unique(N2_euk_bgd$gene) # a vector

GO_terms_euk <- AnnotationDbi::select(GO.db, keys=unique(N2_euk_bgd$GO), columns=c("TERM","DEFINITION"), keytype="GOID") %>% dplyr::rename(TERM = GOID, TERM_NAME = TERM)

# Load in gene set info #

core_test <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/core_tran_WBG.tsv", col_names = c("tran", "gene")) %>% 
  dplyr::select(gene) %>% 
  dplyr::pull() # need to vectorize

acc_test <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/accessory_tran_WBG.tsv", col_names = c("tran", "gene")) %>% 
  dplyr::select(gene) %>%
  dplyr::pull() # need to vectorize

priv_test <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/private_tran_WBG.tsv", col_names = c("tran", "gene")) %>% 
  dplyr::select(gene) %>%
  dplyr::pull() # need to vectorize


# Core Eukaryota - No enrichment...?
enGO_core_euk <- clusterProfiler::enricher(
  gene = core_test,
  TERM2GENE = N2_euk_bgd %>% dplyr::select(GO,gene),
  TERM2NAME = GO_terms_euk %>% dplyr::select(TERM,TERM_NAME),
  universe = N2_euk_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  qvalueCutoff = 0.05,
)

head(enGO_core_euk)

dotplot(enGO_core_euk, showCategory = 40, title = "Core Eukaryota")

# Accessory Eukaryota
enGO_acc_euk <- clusterProfiler::enricher(
  gene = acc_test,
  TERM2GENE = N2_euk_bgd %>% dplyr::select(GO,gene),
  TERM2NAME = GO_terms_euk %>% dplyr::select(TERM,TERM_NAME),
  universe = N2_euk_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  qvalueCutoff = 0.05,
)

head(enGO_acc_euk)

dotplot(enGO_acc_euk, showCategory = 40, title = "Accessory Eukaryota")

# Private Eukaryota - none of the 37 private genes are found in the background set
enGO_priv_euk <- clusterProfiler::enricher(
  gene = priv_test,
  TERM2GENE = N2_euk_bgd %>% dplyr::select(GO,gene),
  TERM2NAME = GO_terms_euk %>% dplyr::select(TERM,TERM_NAME),
  universe = N2_euk_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  qvalueCutoff = 0.05,
)

head(enGO_priv_euk)

dotplot(enGO_priv_euk, showCategory = 40, title = "Private Eukaryota")






### Nematoda ###
N2_nematoda_bgd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output/nematoda_annotations/N2_longestIso_background_nematoda_GOannotations.tsv", col_names = c("transcript","GO")) %>%
  dplyr::left_join(tran_gene, by = "transcript") %>%
  dplyr::select(gene,GO) %>%
  tidyr::separate_rows(GO, sep = ",") 

N2_nematoda_bgd_genes <- unique(N2_nematoda_bgd$gene) # a vector

GO_terms_nematoda <- AnnotationDbi::select(GO.db, keys=unique(N2_nematoda_bgd$GO), columns=c("TERM","DEFINITION"), keytype="GOID") %>% dplyr::rename(TERM = GOID, TERM_NAME = TERM)


# Core Nematoda 
enGO_core_nematoda <- clusterProfiler::enricher(
  gene = core_test,
  TERM2GENE = N2_nematoda_bgd %>% dplyr::select(GO,gene),
  TERM2NAME = GO_terms_nematoda %>% dplyr::select(TERM,TERM_NAME),
  universe = N2_nematoda_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  qvalueCutoff = 0.05,
)

head(enGO_core_nematoda)

dotplot(enGO_core_nematoda, showCategory = 40, title = "Core Nematoda")


# Accessory Nematoda
enGo_acc_nematoda <- clusterProfiler::enricher(
  gene = acc_test,
  TERM2GENE = N2_nematoda_bgd %>% dplyr::select(GO,gene),
  TERM2NAME = GO_terms_nematoda %>% dplyr::select(TERM,TERM_NAME),
  universe = N2_nematoda_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  qvalueCutoff = 0.05,
)

head(enGo_acc_nematoda)

dotplot(enGo_acc_nematoda, showCategory = 40, title = "Accessory Nematoda")


# Private Nematoda - none of the 37 genes in private have GO terms (not in the background set)
enGo_priv_nematoda <- clusterProfiler::enricher(
  gene = priv_test,
  TERM2GENE = N2_nematoda_bgd %>% dplyr::select(GO,gene),
  TERM2NAME = GO_terms_nematoda %>% dplyr::select(TERM,TERM_NAME),
  universe = N2_nematoda_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  qvalueCutoff = 0.05,
)

head(enGo_priv_nematoda)

dotplot(enGo_priv_nematoda, showCategory = 40, title = "Private Nematoda")









