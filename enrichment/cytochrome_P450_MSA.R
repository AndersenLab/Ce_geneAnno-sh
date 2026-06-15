library(dplyr)
library(pheatmap)
library(grid)
library(ggplot2)

# Read single file
pim <- read.table("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/CYTO_P450_SC_OGS/MSA/identity_matrices/OG0003093_identity.tsv", skip = 1) %>%
  dplyr::mutate(V1 = gsub("^[^_]+_([^_]+)_.+$", "\\1",V1))

# Set column names
colnames(pim) <- c("strain", pim$V1)

# Convert to matrix
mat <- as.matrix(pim[, -1])
rownames(mat) <- pim[[1]]
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- pim[[1]]

# Quick check
dim(mat)  # Should be 142 x 142
mat[1:5, 1:5]  # Preview corner

# Plot
pheatmap(
  mat,
  main = "OG0003093",
  color = colorRampPalette(c("gold", "orange", "red"))(100),
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  na_col = "black"
)



# Set paths
matrix_dir <- "/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/CYTO_P450_SC_OGS/MSA/identity_matrices"

# Get all identity matrix files
matrix_files <- list.files(matrix_dir, pattern = "_identity.tsv", full.names = TRUE)

plot_list <- list()
summary_df <- data.frame()

for (file in matrix_files) {
  og_name <- gsub("_identity.tsv$", "", basename(file))
  
  message("Processing: ", og_name)
  
  pim <- read.table(file, skip = 1, fill = TRUE) %>%
    dplyr::mutate(V1 = gsub("^[^_]+_([^_]+)_.+$", "\\1", V1))
  
  # Set column names
  colnames(pim) <- c("strain", pim$V1)
  
  # Convert to matrix
  mat <- as.matrix(pim[, -1])
  rownames(mat) <- pim[[1]]
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- pim[[1]]
  
  min_val <- floor(min(mat, na.rm = TRUE))
  max_val <- ceiling(max(mat, na.rm = TRUE))
  
  # Handle case where all values are identical or very similar
  if (min_val == max_val) {
    min_val <- min_val - 1
  }
  
  message("  Range: ", min_val, " - ", max_val)
  
  # Get lower triangle (excludes diagonal and duplicates)
  lower_tri <- mat[lower.tri(mat)]
  
  # Calculate summary stats
  summary_df <- rbind(summary_df, data.frame(
    orthogroup = og_name,
    mean_identity = mean(lower_tri, na.rm = TRUE),
    sd_identity = sd(lower_tri, na.rm = TRUE),
    median_identity = median(lower_tri, na.rm = TRUE),
    min_identity = min(lower_tri, na.rm = TRUE),
    max_identity = max(lower_tri, na.rm = TRUE),
    n_comparisons = length(lower_tri)
  ))
  
  
  
  grid.newpage()
  plot_list[[og_name]] <- pheatmap(
    mat,
    main = og_name,
    color = colorRampPalette(c("gold", "red", "black"))(100),
    clustering_method = "complete",
    breaks = seq(min_val, max_val, length.out = 51),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    na_col = "black")
}

plot_list[['OG0003616']]
grid.newpage()

for (name in names(plot_list)) {
  grid.newpage()
  print(plot_list[[name]])
  readline(prompt = paste("Viewing:", name, "- Press [Enter] for next..."))
}

### How does percent identity correlate with presence in HDRs across wild strains?
hdr_freq <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/CYTO_P450_SC_OGS/cytoP450_sc_OGs.tsv", col_names = c("orthogroup","HDR_count","found_ever")) %>%
  dplyr::select(-found_ever)

summary_df_hdr <- summary_df %>% dplyr::left_join(hdr_freq, by = "orthogroup")

ggplot(summary_df_hdr, aes(x = reorder(orthogroup, mean_identity), y = mean_identity)) +
  geom_errorbar(aes(ymin = mean_identity - sd_identity, ymax = pmin(mean_identity + sd_identity, 100)), width = 0.2, color = 'orange') +
  geom_point(size = 3, color = 'red') +
  facet_grid(HDR_count ~ ., scales = "free_y", space = "free_y", switch = "y") +
  coord_flip() +
  labs(x = NULL, 
       y = "Mean pairwise identity (%)",
       title = "Sequence identitiy among 142 single-copy orthologs of 38 cytochrome P450 orthogroups",
       subtitle = "Error bars = ± 1 SD") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 16, color = 'black'),
    plot.title = element_text(size = 16, color = 'black'),
    plot.subtitle = element_text(size = 14, color = 'black')) +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.001)))







