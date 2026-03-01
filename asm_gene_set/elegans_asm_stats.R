library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

bad_strains <- c("ECA396","ECA1885", "ECA1887", "ECA2949", "ECA2888", "ECA2291")

stats <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/20251126_master_asmSheet.tsv") %>%
  dplyr::filter(species == "CE") %>%
  dplyr::select(strain,n_contigs,contig_bp,ctg_N50,ctg_L90,ctg_N90,fold_cov,genome_busco) %>%
  dplyr::filter(!strain %in% bad_strains)

plot_df <- stats %>%
  transmute(
    strain,
    contig_bp,
    ctg_N50_mb = ctg_N50 / 1e6,
    ctg_N90_mb = ctg_N90 / 1e6
  ) %>%
  pivot_longer(
    cols = c(contig_bp, ctg_N50_mb, ctg_N90_mb),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("contig_bp", "ctg_N50_mb", "ctg_N90_mb"),
      labels = c("Assembly size (bp)", "Contig N50 (Mb)", "Contig N90 (Mb)")
    )
  )

plot_df <- stats %>% 
  dplyr::transmute(strain,
    `Genome assembly` = contig_bp / 1e6,
    N50 = ctg_N50 / 1e6,
    N90 = ctg_N90 / 1e6) %>%
  tidyr::pivot_longer(cols = c(`Genome assembly`, N50, N90),names_to = "metric",values_to = "value")

top <- ggplot(plot_df, aes(x = metric, y = value)) +
  geom_violin(aes(fill = metric), trim = TRUE, alpha = 0.5) +
  scale_fill_manual(values = c("Genome assembly" = "blue", "N50" = "purple", "N90" = "orange")) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, color="black") +
  geom_point(position = position_jitter(width = 0.12), size = 2.5, color = 'black') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(color = 'gray'),
    panel.border = element_blank(),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.y  = element_text(size = 16, color = 'black'),
    axis.text.x = element_blank(),
    plot.margin = margin(l = 30, t = 2, r =2, b = 30),
    axis.ticks.x  = element_blank()
  ) +
  labs(y = "Megabases") +
  coord_cartesian(ylim = c(100,115))
top


bottom <- ggplot(plot_df, aes(x = metric, y = value)) +
  geom_violin(aes(fill = metric), trim = TRUE, alpha = 0.5) +
  scale_fill_manual(values = c("Genome asembly" = "blue", "N50" = "purple", "N90" = "orange")) +
  stat_summary(fun=mean, geom="crossbar", width=0.5, color="black") +
  geom_point(position = position_jitter(width = 0.12), size = 2.5, color = 'black') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(color = 'gray'),
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    legend.position = 'none',
    plot.margin = margin(l = 30, b = 2, r = 2, t = 30),
    axis.text.x  = element_text(size = 18, color = 'black', face = 'bold'),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 16, color = 'black')
  ) +
  labs(y = "Megabases") +
  coord_cartesian(ylim = c(0,20))
bottom

aligned <- cowplot::align_plots(top, bottom, align = "v", axis = "lr")

base <- cowplot::plot_grid(cowplot::plot_grid(
  aligned[[1]],aligned[[2]],
  nrow = 2) + draw_label("Megabases", x=0, y=0.5, vjust= 1.5, angle=90, size = 18, color = 'black', fontface = 'bold'))

cowplot::ggdraw(base) +
  # Upper break marks
  draw_line(x = c(0.04, 0.06), y = c(0.52, 0.54), size = 1) +
  draw_line(x = c(0.04, 0.06), y = c(0.465, 0.485), size = 1)



