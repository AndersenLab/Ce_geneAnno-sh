library(ggplot2)
library(dplyr)
library(data.table)
library(readr)

pav <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/PAV_CB4856/output/CB4856_over50_PASS_noSNV_variants.tsv", col_names = c("chrom","pos","ref","alt","filter","sv_type","sv_length","strain")) %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::select(chrom, pos, sv_type, sv_length, strain) %>%
  dplyr::mutate(sv_end = pos + sv_length) %>%
  dplyr::mutate(source = "PAV")

syri <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/output/CB4856_SyRI_SVs.tsv") %>%
  dplyr::mutate(sv_length = ifelse(sv_type == "INS", nchar(alt), sv_length)) %>%
  dplyr::filter(sv_length >= 50) %>%
  dplyr::select(chrom, pos, sv_type, sv_length, strain) %>%
  dplyr::mutate(sv_end = pos + sv_length) %>%
  dplyr::mutate(source = "SyRI")

plt_dt <- pav %>% dplyr::bind_rows(syri)

ggplot(plt_dt) +
  geom_point(data = plt_dt %>% dplyr::filter(source == "PAV", sv_type == "DEL"), aes(x = pos / 1e6, y = 1, color = source),size = 1) +
  geom_point(data = plt_dt %>% dplyr::filter(source == "SyRI", sv_type == "DEL"), aes(x = pos / 1e6, y = 0.5, color = source), size = 1) +
  facet_wrap(~chrom, scales = "free_x", ncol = 3) +
  scale_color_manual(values = c("SyRI" = "red", "PAV" = "blue")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 16, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.x = element_text(size = 18, color = 'black')) +
  labs(color = "", x = "Genomic position (Mb)", )



ggplot(plt_dt %>% dplyr::filter(sv_type == "DEL")) + 
  geom_boxplot(aes(x = source, y = sv_length, fill = source), outlier.shape = NA) + 
  geom_jitter(aes(y = sv_length, x = source), color = 'black', width = 0.2) +
  scale_fill_manual(values = c("PAV" = "blue", "SyRI" = "red")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 16, color = 'black'),
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 20, color = 'black', hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 18, color = 'black')) +
  labs(y = "Size", title = "DEL") +
  scale_y_log10()

ggplot(plt_dt %>% dplyr::filter(sv_type == "INS")) + 
  geom_boxplot(aes(x = source, y = sv_length, fill = source), outlier.shape = NA) + 
  geom_jitter(aes(y = sv_length, x = source), color = 'black', width = 0.2) +
  scale_fill_manual(values = c("PAV" = "blue", "SyRI" = "red")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 16, color = 'black'),
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 20, color = 'black', hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 18, color = 'black')) +
  labs(y = "Size", title = "INS") +
  scale_y_log10()

ggplot(plt_dt %>% dplyr::filter(sv_type == "INV")) + 
  geom_boxplot(aes(x = source, y = sv_length, fill = source), outlier.shape = NA) + 
  geom_jitter(aes(y = sv_length, x = source), color = 'black', width = 0.2) +
  scale_fill_manual(values = c("PAV" = "blue", "SyRI" = "red")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 16, color = 'black'),
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 20, color = 'black', hjust = 0.5, face = "bold"),
    axis.title.y = element_text(size = 18, color = 'black')) +
  labs(y = "Size", title = "INV") +
  scale_y_log10()
       
