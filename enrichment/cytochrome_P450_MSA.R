library(dplyr)
library(pheatmap)
library(grid)

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











