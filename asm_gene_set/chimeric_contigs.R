library(ggplot2)
library(readr)
library(dplyr)

soi <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/chimeric_contigs/master_strain_contigs.tsv", col_names = c("strain","contig")) %>%
  dplyr::mutate(strain = ifelse(strain == "PB306", "ECA259", strain))
strains <- soi %>% dplyr::distinct(strain) %>% dplyr::pull(strain)

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) 






chimeric_plotting <- function(strain_of_interest) {
  ctgs <- soi %>% dplyr::filter(strain == strain_of_interest) %>% dplyr::pull(contig)
  
  aln <- nucmer %>%
      dplyr::filter(strain == strain_of_interest & contig %in% ctgs)
  
  plot_chim <- ggplot(aln) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
    facet_wrap(~N2_chr, scales = "free") +
    theme(
      # legend.position = 'none',
      axis.text = element_text(size = 14, color = 'black'),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 16, color = 'black', face = 'bold'),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA),
      plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
    # coord_cartesian(xlim = c(8.36, 8.38)) +
    labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)", title = paste0(strain_of_interest, " chimeric alignment"))
  plot_chim
}







plots <- list()
i=1
for (strain in strains) {
  print(paste0(strain, " ", i, "/", length(strains)))
  i=i+1
  plots[[strain]] <- chimeric_plotting(strain)
}

plots


### N2 chrom III region that splits....
aln <- nucmer %>%
  dplyr::filter(strain == "MY1" & contig == "ptg000006l") %>%
  dplyr::filter(N2_chr == "X") %>%
  dplyr::filter(L1 > 100000)

xmin = min(aln$N2S) / 1e6
xmax = max(aln$N2E)/1e6

plot_chim <- ggplot(aln) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_vline(xintercept = min(aln$N2S) / 1e6) +
  geom_vline(xintercept = max(aln$N2E)/1e6) +
  annotate("text", x = xmin + 0.12, y = 9, label = paste0("min = ", round(xmin, 3), " Mb"), vjust = -0.5, size = 6) +
  annotate("text", x = xmax - 0.12, y = 9, label = paste0("max = ", round(xmax, 3), " Mb"), vjust = -0.5, size = 6) +
  facet_wrap(~N2_chr, scales = "free") +
  theme(
    # legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  # coord_cartesian(xlim = c(8.36, 8.38)) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)", title = "MY1 chimeric alignment")
plot_chim


aln2 <- nucmer %>%
  dplyr::filter(strain == "JU346" & contig == "ptg000001l") %>%
  dplyr::filter(N2_chr == "II") %>%
  dplyr::filter(L1 > 100000)

xmin = min(aln2$N2S) / 1e6
xmax = max(aln2$N2E)/1e6

plot_chim2 <- ggplot(aln2) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_vline(xintercept = min(aln2$N2S) / 1e6) +
  geom_vline(xintercept = max(aln2$N2E)/1e6) +
  annotate("text", x = xmin + 0.12, y = 10, label = paste0("min = ", round(xmin, 3), " Mb"), vjust = -0.5, size = 6) +
  annotate("text", x = xmax - 0.12, y = 10, label = paste0("max = ", round(xmax, 3), " Mb"), vjust = -0.5, size = 6) +
  facet_wrap(~N2_chr, scales = "free") +
  theme(
    # legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  # coord_cartesian(xlim = c(8.36, 8.38)) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)", title = "JU346 chimeric alignment")
plot_chim2



