library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)

# for file in *.vcf; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN' $file | awk -v strain=${file%%.*} -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,strain}'; done 

allcalls <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans/strain_dirs/all_vcfs/141strains_all_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "sv_type","sv_length","strain" ))


filt_calls <- allcalls %>% 
  dplyr::filter(abs(sv_length) >= 50) %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::group_by(sv_type, strain) %>%
  dplyr::mutate(type_count = n()) %>%
  dplyr::ungroup()

levels <- summary <- filt_calls %>%
  dplyr::select(sv_type,strain,type_count) %>%
  dplyr::distinct() %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(desc(type_count)) %>%
  dplyr::ungroup() %>% 
  dplyr::distinct(strain) %>%
  dplyr::pull()

summary <- filt_calls %>%
  dplyr::select(sv_type,strain,type_count) %>%
  dplyr::distinct() %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(desc(type_count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS","DEL","INV"))) %>%
  dplyr::mutate(strain = factor(strain, levels = levels))

count_stats <- summary %>%
  dplyr::group_by(sv_type) %>%
  dplyr::summarise(average_count = mean(type_count)) %>%
  dplyr::mutate(label = paste0(sv_type, " mean: ", round(average_count, 0))) %>%
  dplyr::mutate(ypos = max(summary$type_count) * c(0.95, 0.91, 0.87))

ggplot() +
  geom_col(data = summary, aes(x = strain, y = type_count, fill = sv_type), position = position_dodge(width = 0.8)) +
  geom_hline(data = count_stats, aes(yintercept = average_count, color = sv_type), linetype = "dashed", linewidth = 1) +
  geom_text(data = count_stats, aes(x = Inf, y = ypos, label = label, color = sv_type), hjust = 1.1, vjust = -0.3, size = 6) +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  scale_color_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(y = "Count of SV type for wild strains", fill = "SV type") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 9, angle = 60, hjust = 1),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))










size <- filt_calls %>%
  dplyr::select(sv_type,sv_length) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::arrange(desc(sv_length)) %>%
  dplyr::slice_head()

stats <- filt_calls %>%
  dplyr::select(sv_type,sv_length) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::mutate(average = mean(sv_length)) %>%
  dplyr::distinct(average) %>%
  dplyr::mutate(xpos = ifelse(sv_type == "INV", average + 70000,
                              ifelse(sv_type == "DEL", 27000, 300000))) %>%
  dplyr::mutate(label = paste0(sv_type, " mean: ", round(average, 0))) %>%
  dplyr::mutate(ypos = ifelse(sv_type == "INS", 10000, 
                              ifelse(sv_type == "INV", 150, 10000)))


ggplot() +
  geom_histogram(data = size, aes(x = sv_length, fill = sv_type), binwidth = 300, position = position_dodge(width = 0.8)) +
  geom_vline(data = stats,aes(xintercept = average), color = 'black', linetype = "dashed", linewidth = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  geom_text(data = stats, aes(x = xpos, y = ypos, label = label), color = 'black', hjust = 1.1, vjust = -0.3, size = 6) +
  # scale_color_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  facet_wrap(~sv_type, ncol = 1, scales = "free") +
  labs(y = "Count", x = "SV size (bp)", fill = "SV type") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        legend.position = 'none',
        strip.text.x = element_text(color = 'black', size = 12),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 14),
        axis.title = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legensd for the horizontal lines
  scale_y_log10(expand = c(0,0)) 
  









# Visualizing some of the largest SVs
largest <- filt_calls %>%
  dplyr::select(chrom,pos,sv_type, sv_length, strain) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::arrange(desc(sv_length)) %>%
  dplyr::slice_head(n = 3)

soi <- largest %>% 
  dplyr::distinct(strain) %>%
  dplyr::pull()

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain != "ECA396") %>%
  dplyr::filter(strain %in% soi)


del1 <- nucmer %>% dplyr::filter(strain == "JU2526") %>% dplyr::filter(N2_chr == "V") %>% dplyr::filter(contig == "ptg000001l")
del1_rect <- largest %>% dplyr::filter(strain == "JU2526")

d1 <- ggplot(del1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = del1_rect, aes(xmin = (pos - sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(17.1, 17.5)) +
  ggtitle("Largest deletion in JU2526") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
d1




ins <- nucmer %>% dplyr::filter(strain == "CX11264") %>% dplyr::filter(N2_chr == "V") #%>% dplyr::filter(contig == "ptg000006l")
ins_rect <- largest %>% dplyr::filter(strain == "CX11264")

i1 <- ggplot(ins) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = ins_rect, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.5) +
  theme_bw() +
  theme(
    # legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  # coord_cartesian(xlim = c(19, 22)) +
  ggtitle("Largest insertion in CX11264") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
i1



inv <- nucmer %>% dplyr::filter(strain == "ECA1413" | strain == "ECA2968" | strain == "XZ1515") %>% dplyr::filter(N2_chr == "IV") #%>% dplyr::filter(contig == "ptg000006l")
inv_rect <- largest %>% dplyr::filter(sv_type == "INV")

inv1 <- ggplot(inv) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = inv_rect, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.5) +
  theme_bw() +
  facet_wrap(~strain) +
  theme(
    # legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(6, 8)) +
  ggtitle("Largest inversion") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
inv1


eca <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain == "ECA3088", N2_chr == "V", contig == "ptg000001l")

eca_call <- filt_calls %>% dplyr::filter(strain == "ECA3088", pos > 18000000 & pos < 20000000)

eca1 <- ggplot(eca) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = eca_call %>% dplyr::filter(sv_type != "DEL"), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf, fill = sv_type), alpha = 0.5) +
  geom_rect(data = eca_call %>% dplyr::filter(sv_type == "DEL"), aes(xmin = (pos - sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
  scale_fill_manual(values = c("INS" = "blue", "INV" = "gold")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(16, 21.5), ylim = c(0,7)) +
  ggtitle("ECA3088") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
eca1



# hmm <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
#   dplyr::select(-IDY) %>% dplyr::filter(N2_chr == "IV")
# 
# 
# hmmplot <- ggplot(hmm %>% dplyr::filter(strain == "JU258")) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   # geom_rect(data = eca_call %>% dplyr::filter(sv_type != "DEL"), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf, fill = sv_type), alpha = 0.5) +
#   # scale_fill_manual(values = c("INS" = "blue", "INV" = "gold")) +
#   facet_wrap(~strain, scales = "free")+
#   theme_bw() +
#   theme(
#     # legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA),
#     plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
#   coord_cartesian(xlim = c(15, 19)) +
#   # ggtitle("ECA3088") +
#   labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
# hmmplot
# 
# 
