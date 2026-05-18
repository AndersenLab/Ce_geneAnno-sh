library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)


ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/withN2_BRAKER_orthofinder_run/64_core/OrthoFinder/Results_May16/Orthogroups/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE", N2)) %>%
  dplyr::select(Orthogroup, N2, N2_BRAKER, everything())

strainCol_c2 <- colnames(ortho_genes_dd)

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


private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/withN2_BRAKER_orthofinder_run/64_core/OrthoFinder/Results_May16/Orthogroups/Orthogroups_UnassignedGenes.tsv") %>%
  dplyr::filter(!grepl("MTCE", N2)) %>%
  dplyr::select(Orthogroup, N2, N2_BRAKER, everything())

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


#######################################################################################################################################
### OG gene set classification
#######################################################################################################################################
private_freq = (1 / (length(strainCol_c2_u)))

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
    "private" = "magenta3"), limits = c("core", "accessory", "private"), guide = guide_legend(title = NULL) ) +
  ylab("Orthogroups") + 
  xlab("Genomes") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(classification$sum), by = 25)) +
  theme(
    axis.title = element_text(size = 24, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    plot.margin = margin(l = 20, r = 20, t = 20),
    legend.text = element_text(size=22, color = 'black'),
    axis.text = element_text(size=18, color = 'black'))
gs_allOrtho

OG_class_count <- classification %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_OG = sum(n)) %>%
  dplyr::ungroup()



#######################################################################################################################################
### Plotting horizontal gene set proportion plot
#######################################################################################################################################
classification_genes <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) 
  
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
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  labs(x = "Percent of genes", fill = "Gene set") +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    axis.text.y = element_text(size = 7, color = 'black'),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 19, color = 'black'),
    axis.title.x = element_text(size = 22, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    legend.text = element_text(color = 'black', size = 14),
    legend.title = element_text(color = 'black', size = 16),
    panel.border = element_rect(color = 'black', fill = NA)
  )
contrib


#######################################################################################################################################
### Plotting the number of private genes per strain
#######################################################################################################################################
classification_private <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::filter(class == "private")

private_genes_per_strain <- classification_private %>%
  dplyr::select(-c(sum, freq, class)) %>%  
  tidyr::pivot_longer(cols = -1, 
                      names_to = "strain", 
                      values_to = "presence") %>%
  dplyr::filter(presence == 1) %>%  
  dplyr::count(strain, name = "n_private_genes") %>%
  dplyr::mutate(strain = gsub("_count","",strain)) %>%
  dplyr::arrange(desc(n_private_genes)) %>%
  dplyr::mutate(idy = ifelse(strain == "N2_BRAKER", TRUE,FALSE))


ggplot(private_genes_per_strain, aes(x = reorder(strain, -n_private_genes), y = n_private_genes)) +
  geom_bar(aes(fill = idy), stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "blue ", "FALSE" = "magenta3")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'none',
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.y = element_text(size = 18, color ='black'),
    axis.text.x = element_text(angle = 60, hjust = 1, color = 'black', size = 12)) +
  labs(y = "Number of private genes") +
  scale_y_continuous(expand = c(0,0))
  
  
#######################################################################################################################################
### Number of CNVs between N2 WormBase and N2 BRAKER
#######################################################################################################################################
n2_n2 <- all_relations %>% dplyr::select(N2 = N2_count, N2_BRAKER = N2_BRAKER_count) %>% dplyr::mutate(both_na = ifelse(is.na(N2) & is.na(N2_BRAKER), TRUE, FALSE)) %>%
    dplyr::filter(both_na == FALSE) %>% dplyr::select(-both_na)
  
cnv <- n2_n2 %>% dplyr::filter(!is.na(N2) & !is.na(N2_BRAKER)) %>%
  dplyr::mutate(N2_more = ifelse(N2 > N2_BRAKER,1,0),
                N2_BRAKER_more = ifelse(N2_BRAKER > N2, 1, 0)) %>%
  dplyr::select(N2_more, N2_BRAKER_more) %>%
  dplyr::summarize(across(everything(), ~ sum(.)))
  
stable_copy <- n2_n2 %>% dplyr::filter(!is.na(N2) & !is.na(N2_BRAKER)) %>%
  dplyr::filter(N2 == N2_BRAKER)
  
pav <- n2_n2 %>% dplyr::filter(is.na(N2) | is.na(N2_BRAKER)) %>% dplyr::summarise(across(everything(), ~ sum(!is.na(.))))

pav_N2_WB_ogs <- all_relations %>% dplyr::select(Orthogroup, N2 = N2_count, N2_BRAKER = N2_BRAKER_count) %>% dplyr::mutate(both_na = ifelse(is.na(N2) & is.na(N2_BRAKER), TRUE, FALSE)) %>%
  dplyr::filter(both_na == FALSE) %>% dplyr::select(-both_na) %>%
  dplyr::filter(!is.na(N2) & is.na(N2_BRAKER)) %>% dplyr::pull(Orthogroup)
  
N2_WB_specific_genes <- ortho_genes_dd %>% dplyr::select(Orthogroup, N2) %>% dplyr::filter(Orthogroup %in% pav_N2_WB_ogs)
N2_WB_specific_genes_private <- private_OGs %>% dplyr::select(Orthogroup, N2) %>% dplyr::filter(Orthogroup %in% pav_N2_WB_ogs)

final_N2_WB_specific <- N2_WB_specific_genes %>% dplyr::bind_rows(N2_WB_specific_genes_private) %>% dplyr::select(-Orthogroup) %>%
  tidyr::separate_rows(N2, sep = ", ") %>%
  dplyr::mutate(N2 = gsub("transcript_","", N2)) # 1,031 genes

# write.table(final_N2_WB_specific, "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/withN2_BRAKER_orthofinder_run/N2_WormBase_specific_noN2BRAKER.tsv", row.names = F, quote = F, col.names = F)

# Stats on the 1,031 genes:
wb_1031 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/withN2_BRAKER_orthofinder_run/N2_WormBase_specific_noN2BRAKER_gff.tsv", col_names = c("chrom","source","type","start","end","score","strand","phase","attribute")) %>%
  dplyr::select(-c(source,score,strand,phase)) %>%
  tidyr::separate(attribute, into = c("id","rest"), sep = ";") %>%
  dplyr::mutate(gene = ifelse(type == "mRNA", gsub("ID=transcript:","",id), 
                              ifelse(type == "CDS", gsub("Parent=transcript:","",rest),NA))) %>%
  dplyr::select(-id, -rest)

wb_1031_cds <- wb_1031 %>% dplyr::filter(type == "CDS") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(cds_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene, cds_count) %>%
  dplyr::group_by(cds_count) %>%
  dplyr::mutate(num_genes_cds = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cds_count, num_genes_cds) %>%
  dplyr::mutate(source = "N2 WB specific (not in N2 BRAKER)")

all_wb <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.CDS.nuclear.cleaned.gff3", col_names = c("chrom","source","type","start","end","score","strand","phase","attribute")) %>%
  dplyr::select(-c(source,score,strand,phase,X10)) %>%
  dplyr::group_by(attribute) %>%
  dplyr::mutate(cds_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(attribute,cds_count) %>%
  dplyr::group_by(cds_count) %>%
  dplyr::mutate(num_genes_cds = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(cds_count, num_genes_cds) %>%
  dplyr::mutate(source = 'All N2 genes') %>%
  dplyr::bind_rows(wb_1031_cds)

top <- ggplot(all_wb) +
  geom_col(aes(x = cds_count, y = num_genes_cds, fill = source), position = "dodge") +
  scale_fill_manual(values = c("All N2 genes" = "blue", "N2 WB specific (not in N2 BRAKER)" = "black")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    legend.position = 'inside',
    legend.text = element_text(size = 18, color = 'black'),
    legend.position.inside = c(0.8,0.8)
  ) +
  coord_cartesian(ylim = c(550,3100)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "CDS count per gene", y = "Number of genes", fill = "")
top

bottom <- ggplot(all_wb) +
  geom_col(aes(x = cds_count, y = num_genes_cds, fill = source), position = "dodge") +
  scale_fill_manual(values = c("All N2 genes" = "blue", "N2 WB specific (not in N2 BRAKER)" = "black")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'),
    legend.position = 'none',
    axis.title.x = element_text(size = 16, color = 'black')
  ) +  coord_cartesian(ylim = c(0,430)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "CDS count per gene", y = "Number of genes", fill = "")
bottom

cowplot::plot_grid(
  top, bottom,
  align = "v",
  nrow = 2)


all_wb_sizes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.CDS.nuclear.cleaned.gff3", col_names = c("chrom","source","type","start","end","score","strand","phase","attribute")) %>%
  dplyr::select(-c(source,score,strand,phase,X10)) %>%
  dplyr::group_by(attribute) %>%
  dplyr::mutate(cds_span = end - start) %>%
  dplyr::group_by(attribute) %>%
  dplyr::mutate(cds_total_span = sum(cds_span)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(attribute, cds_total_span) %>%
  dplyr::rename(gene = attribute) %>%
  dplyr::mutate(source = "All N2 genes")


n2_wb_sizes <- wb_1031 %>% dplyr::filter(type == 'CDS') %>%
  dplyr::mutate(cds_span = end - start) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(cds_total_span = sum(cds_span)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(gene, cds_total_span) %>%
  dplyr::mutate(source = "N2 WB specific (not in N2 BRAKER)") %>%
  dplyr::bind_rows(all_wb_sizes) %>%
  dplyr::group_by(source) %>%
  dplyr::mutate(average_size = mean(cds_total_span)) %>%
  dplyr::ungroup()

average_sizes = n2_wb_sizes %>%
  dplyr::distinct(source, average_size)

ggplot(n2_wb_sizes) + 
  geom_histogram(aes(x = cds_total_span / 1e3, fill = source), bins = 500) +
  scale_fill_manual(values = c("All N2 genes" = "blue", "N2 WB specific (not in N2 BRAKER)" = "black")) +
  geom_text(data = average_sizes, 
            aes(x = 7, y = c(800, 700), label = paste0("Mean for ", source, ": ", round(average_size / 1e3, 2), " kb")), 
            size = 5, hjust = 0) +  
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.ticks = element_blank(),
    legend.position = 'inside',
    legend.text = element_text(size = 18, color = 'black'),
    legend.position.inside = c(0.8,0.8),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 16, color = 'black')) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "CDS span (kb)", y = "Count", fill = "") +
  coord_cartesian(xlim = c(0,10))



# N2 BRAKER specific (not found in N2 WormBase)
pav_N2_BRAKER_ogs <- all_relations %>% dplyr::select(Orthogroup, N2 = N2_count, N2_BRAKER = N2_BRAKER_count) %>% dplyr::mutate(both_na = ifelse(is.na(N2) & is.na(N2_BRAKER), TRUE, FALSE)) %>%
  dplyr::filter(both_na == FALSE) %>% dplyr::select(-both_na) %>%
  dplyr::filter(is.na(N2) & !is.na(N2_BRAKER)) %>% dplyr::pull(Orthogroup)

N2_BRAKER_specific_genes <- ortho_genes_dd %>% dplyr::select(Orthogroup, N2_BRAKER) %>% dplyr::filter(Orthogroup %in% pav_N2_BRAKER_ogs)
N2_BRAKER_specific_genes_private <- private_OGs %>% dplyr::select(Orthogroup, N2_BRAKER) %>% dplyr::filter(Orthogroup %in% pav_N2_BRAKER_ogs)

final_N2_BRAKER_specific <- N2_BRAKER_specific_genes %>% dplyr::bind_rows(N2_BRAKER_specific_genes_private) %>% dplyr::select(-Orthogroup) %>%
  tidyr::separate_rows(N2_BRAKER, sep = ", ") # 3,492 genes



#######################################################################################################################################
### Looking at N2-specific ###
#######################################################################################################################################
all_relations_clean <- all_relations %>% dplyr::rename_with(~ gsub("_count", "", .x))
  
names_all <- colnames(all_relations_clean %>% dplyr::select(-N2, -Orthogroup))
  
nonRefGenes = as.data.frame(matrix(ncol = 3, nrow = 141))
colnames(nonRefGenes) <- c("strain","N2_specific_genes",'N2_specific_Orthogroups')

OG_list <- list()
  
for (i in 1:length(names_all)) {
    soi <- names_all[i]
    print(paste0("On strain: ", soi, ". ", i, "/142."))
    
    nonRefGenes[i,1] = soi
    
    non_ref <- all_relations_clean %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% dplyr::select(N2) %>% sum(na.rm = TRUE) 
    
    nonRefGenes[i,2] = non_ref
    
    non_ref_og_count <- all_relations_clean %>% dplyr::select(N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% nrow()
    non_ref_og <- all_relations_clean %>% dplyr::select(Orthogroup, N2, .data[[soi]]) %>% dplyr::filter(!is.na(N2) & is.na(.data[[soi]])) %>% dplyr::pull(Orthogroup)
    OG_list[[i]] <- non_ref_og
    
    nonRefGenes[i,3] = non_ref_og_count
    
}
  
nonRefGenes_long <- nonRefGenes %>%
  pivot_longer(
    cols = c(N2_specific_genes, N2_specific_Orthogroups),
    names_to = "metric",
    values_to = "count")
  
label_df_top <- nonRefGenes_long %>%
  dplyr::filter(metric == "N2_specific_genes")  %>% 
  dplyr::arrange(desc(count)) %>% dplyr::slice_head(n = 10)
  
label_df_bottom <- nonRefGenes_long %>%
  dplyr::filter(metric == "N2_specific_genes") %>% 
  dplyr::arrange(count) %>% dplyr::slice_head(n = 10)
  
labels_df <- label_df_top %>% dplyr::bind_rows(label_df_bottom)
  
n2_spec <- ggplot(nonRefGenes_long, aes(x = metric, y = count)) +
geom_boxplot(aes(fill = metric), outlier.size = 0.6, width = 0.7, outlier.shape = NA, alpha = 0.5) +
geom_line(aes(group = strain), alpha = 0.3) +
  geom_point(aes(group = strain),size = 1.5, alpha = 0.6) +
  geom_text_repel(data = labels_df, aes(label = strain), size = 4, max.overlaps = 100) +
  labs(y = "Count") +
  scale_fill_manual(values = c("N2_specific_genes" = "orange", "N2_specific_Orthogroups" = "#DB6333")) +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black')) 
n2_spec
  
# Looking at proprotion in each gene set
OG_vector_N2_spec <- unique(unlist(OG_list))

OG_classes_spec <- all_relations %>%
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
  ) %>% dplyr::select(Orthogroup, class) %>%
  dplyr::filter(Orthogroup %in% OG_vector_N2_spec)

OG_propN2_spec <- OG_classes_spec %>%
  dplyr::count(class, name = "class_count") %>%
  dplyr::mutate(prop = class_count / sum(class_count) * 100) %>%
  dplyr::mutate(class = factor(class, levels = c("accessory", "core", "private"))) %>%
  dplyr::mutate(source = "N2")
  
n2_spec_prop <- ggplot(OG_propN2_spec, aes(x = prop, y = "", fill = class)) +
  geom_col(position = "stack", width = 0.5, color = "black") +
  scale_fill_manual(values = c("accessory" = "#DB6333", "private" = "magenta3", "core" = "green4")) +
  labs(x = "Proportion (%)", y = NULL, fill = "Gene set") +
  theme(
    panel.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_text(size = 16, color = 'black')) 

n2_specific_genes <- cowplot::plot_grid(
  n2_spec, n2_spec_prop,
  nrow = 2,
  rel_heights = c(6,1))
n2_specific_genes
  
  
  
  
#######################################################################################################################################
### New OG classification
#######################################################################################################################################
colnames(all_relations) <- strainCol_c2

n_strains <- length(strainCol_c2_u)
n_wild_strains <- n_strains - 2  # excluding N2 WB and N2_BRAKER

new_class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(
    sum = rowSums(across(-1), na.rm = TRUE),
    sum_wild = rowSums(across(-c(1, N2, N2_BRAKER, sum)), na.rm = TRUE)) %>%
  dplyr::mutate(freq = sum / n_strains) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private")) %>%
  dplyr::mutate(
    updated_class = case_when(
      !is.na(N2) & is.na(N2_BRAKER) & sum_wild == 0 ~ "private_N2_WB", # private N2 WormBase (only in N2 WB, not in N2_BRAKER or any wild strain)
      is.na(N2) & !is.na(N2_BRAKER) & sum_wild == 0 ~ "private_N2_BRAKER", # private N2 BRAKER (only in N2_BRAKER, not in N2 WB or any wild strain)
      !is.na(N2) & is.na(N2_BRAKER) & sum_wild == n_wild_strains ~ "core_N2_WB", # core N2 WormBase (found in all except N2_BRAKER)
      is.na(N2) & !is.na(N2_BRAKER) & sum_wild == n_wild_strains ~ "core_N2_BRAKER", # core N2 BRAKER (found in all except N2 WormBase)
      !is.na(N2) & is.na(N2_BRAKER) & sum_wild > 0 & sum_wild < n_wild_strains ~ "accessory_N2_WB", # accessory N2 WormBase (in N2 WB, not N2_BRAKER, in some but not all wild strains)
      is.na(N2) & !is.na(N2_BRAKER) & sum_wild > 0 & sum_wild < n_wild_strains ~ "accessory_N2_BRAKER", # accessory N2 BRAKER (in N2_BRAKER, not N2, in some but not all wild strains)
      is.na(N2) & is.na(N2_BRAKER) & sum_wild == n_wild_strains ~ "core_WSs", # found in all wild strains but not N2 WB or N2 BRAKER
      is.na(N2) & is.na(N2_BRAKER) & sum_wild > 0 & sum_wild < n_wild_strains & class != "private" ~ "accessory_WSs", # found in some wild strains but not N2 WB or N2 BRAKER
      TRUE ~ class # otherwise the OG class stays the same
    )
  )

class_stats <- new_class %>%
  dplyr::select(Orthogroup, N2, N2_BRAKER, sum, sum_wild, freq, class, updated_class) %>% 
  dplyr::group_by(class) %>%
  dplyr::mutate(OG_class_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(updated_class) %>%
  dplyr::mutate(updated_class_count = n()) %>%
  dplyr::ungroup() 

acc <- class_stats %>% dplyr::filter(class == "accessory") %>% dplyr::select(class, updated_class, OG_class_count, updated_class_count) %>% dplyr::distinct()
core <- class_stats %>% dplyr::filter(class == "core") %>% dplyr::select(class, updated_class, OG_class_count, updated_class_count) %>% dplyr::distinct()
private <- class_stats %>% dplyr::filter(class == "private") %>% dplyr::select(class, updated_class, OG_class_count, updated_class_count) %>% dplyr::distinct()

cat <- core %>% 
  dplyr::bind_rows(acc, private) %>% 
  dplyr::mutate(proportion = updated_class_count / OG_class_count) %>%
  dplyr::mutate(updated_class = ifelse(updated_class == "core_WSs", "core_WSs (2)", 
                                       ifelse(updated_class == "core_N2_WB", "core_N2_WB (1)", updated_class))) %>%
  dplyr::arrange(class, desc(proportion)) %>%
  dplyr::mutate(updated_class = factor(updated_class, levels = unique(updated_class))) %>%
  dplyr::mutate(class = factor(class, levels = c("core", "accessory", "private"))) 

class_colors <- c(
  "core" = "green4",
  "accessory" = "#DB6333",
  "accessory_WSs" = "darkorange",
  "accessory_N2_BRAKER" = "chocolate4",
  "accessory_N2_WB" = "tan3",
  "core_N2_BRAKER" = "orange",
  "core_N2_WB (1)" = "firebrick3",
  "core_WSs (2)" = "brown",
  "private" = "magenta3",
  "private_N2_WB" = "violet",
  "private_N2_BRAKER" = "purple")

ggplot(cat, aes(x = class, y = proportion, fill = updated_class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = class_colors) +
  theme(
    axis.text.x = element_text(size = 20, color = 'black'),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.title = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 16, color = 'black')) +
  labs(y = "Proportion", fill = "Updated class")


### Looking at gene counts!
keys <- new_class %>% dplyr::select(Orthogroup, class, updated_class)

all_relations_counts <- all_relations %>% dplyr::left_join(keys, by = "Orthogroup") %>%
  dplyr::mutate(count = rowSums(across(-c(Orthogroup, class, updated_class)), na.rm = TRUE))

updated_class_gene_count <- all_relations_counts %>% dplyr::select(updated_class, count) %>%
  dplyr::group_by(updated_class) %>%
  dplyr::mutate(gene_count = sum(count)) %>%
  dplyr::ungroup() %>% 
  dplyr::distinct(updated_class, gene_count) %>%
  dplyr::mutate(updated_class = ifelse(updated_class == "core_WSs", "core_WSs (2)", 
                                       ifelse(updated_class == "core_N2_WB", "core_N2_WB (1)", updated_class))) 

cat_count <- cat %>% 
  dplyr::left_join(updated_class_gene_count, by = "updated_class") %>%
  dplyr::arrange(class, desc(proportion)) %>%
  dplyr::mutate(updated_class = factor(updated_class, levels = unique(updated_class))) %>%
  dplyr::mutate(class = factor(class, levels = c("core", "accessory", "private"))) 

# Adding gene counts for each updated gene set to the plot
ggplot(cat_count, aes(x = class, y = proportion, fill = updated_class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(label = scales::comma(gene_count)), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = class_colors) +
  theme(
    axis.text.x = element_text(size = 20, color = 'black'),
    axis.text.y = element_text(size = 16, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    legend.title = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 16, color = 'black')) +
  labs(y = "Proportion", fill = "Updated class")


single_copy_allButBRAKER_og <- new_class %>% dplyr::filter(updated_class == 'core_N2_WB') %>% dplyr::pull(Orthogroup)
single_copy_allButBRAKER <- ortho_genes_dd %>% dplyr::filter(Orthogroup %in% single_copy_allButBRAKER) %>% dplyr::pull(N2)



