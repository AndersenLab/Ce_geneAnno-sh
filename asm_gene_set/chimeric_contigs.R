library(ggplot2)
library(readr)
library(dplyr)


# MY1
nuc <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::filter(strain == 'MY1') %>% dplyr::filter(N2_chr == "II" | N2_chr == "III" | N2_chr == "X", contig == "ptg000006l" | contig == "ptg000013l")

conservedINV <- ggplot(nuc) +
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
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
conservedINV