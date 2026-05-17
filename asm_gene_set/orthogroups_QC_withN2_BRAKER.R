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
  dplyr::filter(!grepl("MTCE",N2))

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
  geom_bar(stat = "identity", alpha = 0.5, color = "black", linewidth = 0.3) +
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
  new_class <- all_relations %>%
    dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
    dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
    dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
    dplyr::mutate(
      class = case_when(
        freq == 1 ~ "core",
        freq > private_freq & freq < 1 ~ "accessory",
        freq == private_freq ~ "private",
        TRUE ~ "undefined")) 
  
  
