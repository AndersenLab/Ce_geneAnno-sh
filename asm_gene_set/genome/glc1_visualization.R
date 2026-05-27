library(dplyr)
library(readr)
library(ggplot2)

df = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/elegans.merged.1kbCOV.5kbALIGN.annotated.final.vcf")

CB4856var <- df %>%
  dplyr::rename(CHROM = '#CHROM', start = POS) %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856) %>%
  dplyr::filter(CHROM == "V", (start >= 16115967 & start <= 16276907)) %>%
  dplyr::filter(CB4856 != "./.") 

plot_df <- CB4856var %>%
  dplyr::mutate(
    lenDEL = case_when(INFO == "DEL" ~ (end - start), TRUE ~ NA_real_),
    lenINS = case_when(INFO == "INS" ~ (end - start), TRUE ~ NA_real_))

SNPs <- plot_df %>%
  dplyr::filter(INFO == "SNP") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856) %>%
  dplyr::filter(!grepl(",", ALT))

deletions <- plot_df %>%
  dplyr::filter(INFO == "DEL") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenDEL) %>%
  dplyr::filter(lenDEL >= 50)
  
insertions <- plot_df %>%
  dplyr::filter(INFO == "INS") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenINS) %>%
  dplyr::filter(lenINS >= 50)

# Load in SNPs called by GATK pipeline with SR data
SR = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/SR.GATK.CB4856.final.vcf")

GATK <- SR %>%
  dplyr::filter(POS >= 16115967 & POS <= 16276907) %>%
  dplyr::rename(start = POS) %>%
  dplyr::mutate(end = start+1) %>%
  dplyr::filter(CB4856 != './.' & CB4856 != '0/0')

indels = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/SR.GATK.CB4856.indels.final.vcf")

G_del <- indels %>%
  dplyr::filter(START >= 16115967 & START <= 16276907) %>%
  dplyr::filter(ANNO == "DEL") %>%
  dplyr::filter() %>%
  dplyr::filter((END - START) < 50)
  
G_ins <- indels %>%
  dplyr::filter(START >= 16115967 & START <= 16276907) %>%
  dplyr::filter(ANNO == "INS") %>%
  dplyr::filter((END - START) < 50)




hist_data <- ggplot_build(
  ggplot(GATK, aes(x = start / 1e6)) + geom_histogram(binwidth = 0.0002)
)$data[[1]]

y_max <- max(hist_data$count)
print(y_max)
y_scaler <- 0.2 / y_max

hist <- ggplot(GATK, aes(x = start / 1e6)) +
  geom_histogram(binwidth = 0.0002, aes(y = ..count.. * y_scaler, fill = 'SNVs (GATK)'), alpha = 0.8) + 
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  scale_fill_manual(
    name = "",
    values = c("SNVs (GATK)" = "purple"),
    breaks = c("SNPs (GATK)")
  ) +
  scale_y_continuous(limits = c(0, 0.2)) 

hist



plot <- ggplot() +
  geom_rect(data = plot_df, aes(xmin = min(start)/1e6, xmax = max(start)/1e6, ymin = 1.75, ymax = 2.25, fill = 'QTL')) +
  geom_rect(data = plot_df, aes(xmin = 16219634/1e6, xmax = 16221917/1e6, ymin = 1.75, ymax = 2.25, fill = 'glc-1')) +

  geom_rect(data = hist_data, aes(xmin = xmin - 0.0001, xmax = xmin + 0.0001, ymin = 1.5,  ymax = count * y_scaler + 1.5, fill = 'SNVs (GATK)'), alpha = 0.8) +
  
  geom_rect(data = GATK, aes(xmin = start/1e6 - 0.00001, xmax = end/1e6, ymin = 1.0, ymax = 1.5, fill = 'SNVs (GATK)')) + 
  geom_rect(data = deletions, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.25, ymax = 0.75, fill = 'Deletions')) +
  geom_rect(data = G_del, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.25, ymax = 0.75, fill = 'Deletions')) +
  geom_rect(data = G_ins, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.25, ymax = 0.75, fill = 'Insertions')) +
  geom_rect(data = insertions, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.25, ymax = 0.75, fill = 'Insertions')) +
  
  ggtitle("V") +
  
  scale_x_continuous(name = "Genomic Position (Mb)", labels = scales::number_format(scale = 1, accuracy = 0.01)) +
  
  scale_fill_manual(
    name = "",
    values = c("QTL" = "darkolivegreen", 'glc-1' = 'darkolivegreen3', "SNVs (GATK)" = "purple", "Insertions" = "blue", "Deletions" = "red"),
    breaks = c("QTL", "glc-1", "SNPs (GATK)", "Insertions", "Deletions")
  ) + 
  
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.background = element_blank()
  ) 
plot



ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc_1.jpg", plot, dpi=900, width = 7.5, height = 6)













del_all <- plot_df %>%
  dplyr::filter(INFO == "DEL") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenDEL) 
ins_all <- plot_df %>%
  dplyr::filter(INFO == "INS") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenINS) 



glc <- ggplot() +
  # Plotting common, narrowed QTL
  geom_rect(data = plot_df, aes(xmin = 16.21, xmax = 16.23, ymin = 3.75, ymax = 4.25, fill = 'QTL')) +
  geom_rect(data = plot_df, aes(xmin = 16219634/1e6, xmax = 16221917/1e6, ymin = 3.75, ymax = 4.25, fill = 'glc-1')) +
  # Plotting variants
  geom_rect(data = SNPs, aes(xmin = start/1e6 - 0.0001, xmax = end/1e6 + 0.0001, ymin = 2.75, ymax = 3.25, fill = 'SNPs (paftools)')) + 
  geom_point(data = SNPs, aes(x = (start + end) / 2 / 1e6, y = 2.75), size = 0.5) +
  geom_rect(data = GATK, aes(xmin = start/1e6 - 0.0001, xmax = end/1e6 + 0.0001, ymin = 1.75, ymax = 2.25, fill = 'SNPs (GATK)')) + 
  geom_point(data = GATK, aes(x = (start + end) / 2 / 1e6, y = 1.75), size = 0.5) +
  geom_rect(data = del_all, aes(xmin = start/1e6 - 0.0001, xmax = end/1e6 + 0.0001, ymin = 0.75, ymax = 1.25, fill = 'Deletions')) +
  geom_point(data = del_all, aes(x = (start + end) / 2 / 1e6, y = 1.25), size = 0.5) +
  geom_rect(data = ins_all, aes(xmin = start/1e6 - 0.00001, xmax = end/1e6 + 0.00001, ymin = 0.75, ymax = 1.25, fill = 'Insertions')) +
  geom_point(data = ins_all, aes(x = (start + end) / 2 / 1e6, y = 0.75), size = 0.5) +
  
  
  ggtitle("V") +
  
  scale_x_continuous(name = "Genomic Position (Mb)", labels = scales::number_format(scale = 1, accuracy = 0.01), limits = c(16.21, 16.23)) +
  
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black"),  
    axis.title.x = element_text(size = 14, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5), 
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.background = element_rect(color = "white", fill = "white")
  ) +
  
  
  scale_fill_manual(
    name = "", 
    values = c("QTL" = "darkolivegreen", 'glc-1' = 'darkolivegreen3', "SNPs (paftools)" = "purple", "SNPs (GATK)" = "darkorange", "Insertions" = "blue", "Deletions" = "red"),
    breaks = c("QTL", "glc-1", "SNPs (paftools)", "SNPs (GATK)", "Insertions", "Deletions")
  )

glc




total_gatk <- nrow(GATK)
total_paf <- nrow(SNPs) # right now paftools has muli-allelic sites, I need to filter to bi-allelic
print(nrow(deletions))
print(nrow(insertions))

# Shared SNV calls between paftools and GATK
shared <- SNPs %>%
  dplyr::filter(start %in% GATK$start)

# Unique to paftools
unique_SNPs <- SNPs %>%
  dplyr::filter(!(start %in% GATK$start))

# Unique to GATK 
unique_GATK <- GATK %>%
  dplyr::filter(!(start %in% SNPs$start))

num_shared <- nrow(shared)
num_unique_SNPs <- nrow(unique_SNPs)
num_unique_GATK <- nrow(unique_GATK)

cat("Number of shared SNVs:", num_shared, "\n")
cat("Unique to paftools:", num_unique_SNPs, "Total:", total_paf, "\n")
cat("Unique to GATK:", num_unique_GATK, "Total:", total_gatk, "\n")