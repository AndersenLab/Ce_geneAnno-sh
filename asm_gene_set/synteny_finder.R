library(ape)
library(dplyr)
library(ggplot2)

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
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=16, color = 'black'),
    axis.text = element_text(size=14, color = 'black')
  )
gs_allOrtho

HOG_class_count <- classification %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_HOG = sum(n)) %>%
  dplyr::ungroup()


# Splitting complex orthologous relationships into syntenic 1to1 orthologs
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
all_genes_strain <- genes_strain %>%
  dplyr::filter(strain != "N2" | grepl("protein_coding", attributes)) %>% #filter out non-protein coding N2 genes
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) 

# test1 <- all_genes_strain %>%
#   dplyr::filter(strain == "N2") # ~19,000

# write.table(all_genes_strain,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/116strain_genes.tsv", quote = F, row.names = F, col.names = T, sep = '\t')


nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/nucmer_runs/115_WI_transformed.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-L1,-L2,-IDY,-LENR,-LENQ)
write.table(nucmer,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/115WS_nucmer_clean.tsv", quote = F, row.names = F, col.names = T, sep = '\t')