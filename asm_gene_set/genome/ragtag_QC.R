library(ggplot2)
library(readr)
library(dplyr)


asm20 <- read.table("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/asm20_ragtag/unplaced_ctg_seq.tsv", sep = ' ', col.names = c("strain","asm20_number","asm20_seq")) %>%
  dplyr::filter(strain != "ECA2555")
 
asm5 <- read.table("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/asm5_default/scaffolded_asm/asm5_ragtag/unplaced_ctg_seq.tsv", sep = ' ', col.names = c("strain","asm5_number","asm5_seq"))

cat <- asm20 %>% dplyr::left_join(asm5, by = 'strain')


df_long <- cat %>%
  tidyr::pivot_longer(
    cols = c(asm20_seq, asm5_seq, asm20_number, asm5_number),
    names_to = c("divergence", ".value"),
    names_pattern = "(asm[0-9]+)_(.*)")

strain_order <- cat %>%
  dplyr::arrange(asm5_seq) %>%
  dplyr::pull(strain)

df_long$strain <- factor(df_long$strain, levels = strain_order)

segments_df <- df_long %>%
  dplyr::select(strain, divergence, seq) %>%
  tidyr::pivot_wider(names_from = divergence, values_from = seq)

ggplot(df_long, aes(x = strain, y = seq / 1e6, color = divergence, size = number)) +
  geom_segment(data = segments_df, aes(x = strain, xend = strain, y = asm5 / 1e6, yend = asm20 / 1e6), linetype = "dotted", color = "gray40", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point() +
  scale_color_manual(values = c("asm20" = "red", "asm5" = "blue")) +
  theme_bw() +
  labs(
    x = "Strain",
    y = "Unplaced sequence (Mb)",
    size = "Number of unplaced contigs",
    color = "Assembly method",
    title = "Scaffolding using CGC1"
  ) +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, color = 'black'),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    legend.position = 'inside',
    plot.title = element_text(size = 15, color= 'black', hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = 'black'),
    legend.position.inside = c(0.6,0.8),
    legend.text = element_text(size = 12, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 16, color = 'black'))

