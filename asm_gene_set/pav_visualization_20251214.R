library(readr)
library(plyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(circlize)
library(data.table)
library(ComplexHeatmap)  
library(grid)           

# for file in *.vcf; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN' $file | awk -v strain=${file%%.*} -v OFS='\t' '$7 >= 50 || $7 <= -50 {print $1,$2,$3,$4,$5,$6,$7,strain}'; done | grep -w "PASS" | grep -v -w "SNV"

# allcalls <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans/strain_dirs/all_vcfs/141strains_all_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "sv_type","sv_length","strain" )) # did not include the 'PASS' filter - contains compound variants (variants contianed insidce larger variants)
allcalls <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/elegans/strain_dirs/all_vcfs/141_over50_PASS_variants.tsv", col_names = c("chrom", "pos", "ref", "alt", "filter", "sv_type","sv_length","strain")) %>% dplyr::select(-filter)

filt_calls <- allcalls %>% 
  # dplyr::filter(abs(sv_length) >= 50) %>%
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




prop <- summary %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(total = sum(type_count)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sv_type,strain) %>%
  dplyr::mutate(prop = type_count / total * 100) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INV","DEL","INS"))) 

ggplot() +
  geom_bar(data = prop, aes(x = strain, y = prop, fill = sv_type), stat = "identity") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(y = "Proportion of SVs", fill = "SV type") +
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

levels2 <- filt_calls %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(total_SV_length = sum(sv_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,total_SV_length) %>%
  dplyr::arrange(desc(total_SV_length)) %>%
  dplyr::distinct(strain) %>%
  dplyr::pull()

total_size <- filt_calls %>%
  dplyr::select(sv_type,sv_length,strain) %>%
  dplyr::group_by(strain,sv_type) %>%
  dplyr::mutate(type_total_length = sum(sv_length)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,sv_type, type_total_length) %>%
  dplyr::mutate(strain = factor(strain, levels = levels2)) %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INV","DEL","INS"))) #%>%
  # dplyr::group_by(strain) %>%
  # dplyr::mutate(total_SV_length = sum(type_total_length)) %>%
  # dplyr::ungroup() %>%
  # dplyr::arrange(desc(total_SV_length)) %>%
  # dplyr::filter(strain == "ECA3088" | strain == "ECA2199" | strain == "ECA1493")


ggplot() +
  geom_bar(data = total_size, aes(x = strain, y = type_total_length / 1e6, fill = sv_type), stat = "identity") +
  scale_fill_manual(values = c("INS" = "blue", "DEL" = "red", "INV" = "gold")) +
  labs(y = "Length of SVs (Mb)", fill = "SV type") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 9, angle = 60, hjust = 1),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  coord_cartesian(ylim = c(0,14)) +
  scale_y_continuous(breaks = seq(0,14, by =2), expand = c(0,0))







size <- filt_calls %>%
  dplyr::select(sv_type,sv_length) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::arrange(desc(sv_length))# %>%
  # dplyr::slice_head()

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
  dplyr::slice_head(n = 4) #%>%
  dplyr::filter(sv_length == "108633")

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
  geom_rect(data = del1_rect, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
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
  coord_cartesian(xlim = c(17.1, 18)) +
  ggtitle("Largest deletion in JU2526 (V)") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
d1



del2 <- ggplot(nucmer %>% dplyr::filter(N2_chr == "V")) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = largest, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
  theme_bw() +
  facet_wrap(~strain) +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(2.1, 2.9)) +
  ggtitle("Second largest deletion (V)") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
del2





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
inv_rect <- largest %>% dplyr::filter(sv_type == "INV") %>% dplyr::filter(strain != "JU3226")

inv1 <- ggplot(inv) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = inv_rect, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
  theme_bw() +
  facet_wrap(~strain) +
  theme(
    legend.position = 'none',
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


N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") 
n2_genes <- N2_gff %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::select(seqid,start,end,attributes)

inv1gm <- ggplot(inv) +
  geom_rect(data = n2_genes %>% dplyr::filter(seqid == "IV"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = "green", alpha = 0.3) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme_bw() +
  facet_wrap(~strain) +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  # coord_cartesian(xlim = c(6.8, 7.5), ylim = c(5,13)) +
  ggtitle("Largest inversion (IV)") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
inv1gm


eca <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain == "ECA3088", N2_chr == "V", contig == "ptg000001l")

eca_INV <- eca %>% dplyr::filter(N2E > N2S)

eca_call <- filt_calls %>% dplyr::filter(strain == "ECA3088", chrom == 'V', pos > 18000000 & pos < 20000000)

eca1 <- ggplot(eca) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  # geom_rect(data = eca_call %>% dplyr::filter(sv_type == "DEL"), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
  # geom_rect(data = eca_call %>% dplyr::filter(sv_type != "DEL"), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf, fill = sv_type), alpha = 0.5) +
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
  coord_cartesian(xlim = c(18.36, 19.70), ylim = c(0,7)) +
  ggtitle("ECA3088") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
eca1


eca_all <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain == "ECA3088")
  
eca_allplot <- ggplot(eca_all) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  scale_fill_manual(values = c("INS" = "blue", "INV" = "gold")) +
  facet_wrap(~N2_chr, scales = 'free') +
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
  # coord_cartesian(xlim = c(17, 21), ylim = c(0,7)) +
  ggtitle("ECA3088") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
eca_allplot


eca_callFour <- filt_calls %>% dplyr::filter(strain == "ECA3088", chrom == "IV", pos > 13000000 & pos < 18000000)

ecaFour <- ggplot(eca_all %>% dplyr::filter(N2_chr == "IV", contig == "ptg000002l")) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = eca_callFour %>% dplyr::filter(sv_type == "DEL"), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
  geom_rect(data = eca_callFour %>% dplyr::filter(sv_type != "DEL"), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf, fill = sv_type), alpha = 0.5) +
  scale_fill_manual(values = c("INS" = "blue", "INV" = "gold")) +
  facet_wrap(~N2_chr, scales = 'free') +
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
  coord_cartesian(xlim = c(13, 18)) +
  ggtitle("ECA3088") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
ecaFour



eca <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(N2_chr == "IV") %>% dplyr::filter(strain != "ECA396")

ecaAll <- ggplot(eca) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  facet_wrap(~strain, scales = 'free') +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(12,18)) +
  # ggtitle("ECA3088") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
ecaAll
# large inversion on right arm of chrosome form might also be in: ECA1997, ECA2151, ECA2199, ECA2473, ECA248, ECA259, JU311. NIC1790, QX1791
# - no WSs have the inversion in the same contig - a lof of contig breaks in the same locus and then an inverted alignment of a new contig


hm <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(N2_chr == "V") %>% dplyr::filter(strain != "ECA396")

hm1 <- ggplot(hm) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  facet_wrap(~strain) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(18, 20), ylim = c(0,7)) +
  ggtitle("CGC1") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
hm1


hmm <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(N2_chr == "V") %>% dplyr::filter(strain != "ECA396")

hmm2 <- ggplot(hmm) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  facet_wrap(~strain) +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(18, 20), ylim = c(0,7)) +
  ggtitle("CGC1") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
hmm2



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




# Looking at an INV that is found among 
threeINV <- filt_calls %>% 
  dplyr::filter(chrom == "IV", sv_type == "INV") %>%
  dplyr::group_by(pos) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(count != "1") %>%
  dplyr::group_by(pos) %>%
  dplyr::mutate(pos_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(pos_count == max(pos_count))

strains <- threeINV %>%
  dplyr::distinct(strain) %>%
  dplyr::pull()

threeINVnucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::filter(strain %in% strains) %>% dplyr::filter(N2_chr == "IV") %>%
  dplyr::filter(N2S < 8373000 & N2E < 8375200 & N2E > 8373000 | 
                  N2S > 8373000 & N2S < 8375200 & N2E > 8375200 | 
                  N2S > 8373000 & N2E < 8375200) %>%
  dplyr::group_by(strain, contig) %>%
  dplyr::mutate(summedL2 = sum(L2)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(strain) %>%
  dplyr::filter(summedL2 == max(summedL2)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(slope = ifelse(N2S < 8373000 | N2E > 8375200, (WSS - WSE)/(N2S-N2E), NA)) %>%
  dplyr::mutate(inv = ifelse(WSE < WSS, T, F)) %>%
  dplyr::mutate(intercept = WSS - (slope * N2S)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(new_start = ifelse(N2S < 8373000 & N2E < 8375200 & inv == F, ((slope * 8373000) + intercept), 
                                   ifelse(N2S > 8373000 & N2E > 8375200  & inv == T, ((slope * 8373000) + intercept), WSS))) %>%
  dplyr::mutate(new_end = ifelse(N2S > 8373000 & N2E > 8375200 & inv == F, ((slope * 8375200) + intercept), 
                                 ifelse(N2S < 8373000 & N2E < 8375200 & inv == T, ((slope * 8373000) + intercept), WSE))) %>%
  dplyr::select(-IDY,-L2,-L1,-LENR,-LENQ)


conservedINV <- ggplot(threeINVnucmer) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = new_start / 1e6, yend = new_end / 1e6, color = contig), linewidth = 1) +
  geom_rect(data = threeINV %>% dplyr::slice_head(n=1) %>% dplyr::select(chrom,pos,sv_length), aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
  facet_wrap(~strain, scales = "free") +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)) +
  coord_cartesian(xlim = c(8.363, 8.3852)) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
conservedINV







win_start <- 8373000
win_end   <- 8375200

threeINVnucmer_trimmed <- threeINVnucmer %>%
  # keep only alignments that overlap the window at all
  dplyr::filter(N2E >= win_start, N2S <= win_end) %>%
  # normalize so N2_start < N2_end, and match W coords to that direction
  dplyr::mutate(
    N2_start = pmin(N2S, N2E),
    N2_end   = pmax(N2S, N2E),
    W_start  = if_else(N2S <= N2E, WSS, WSE),
    W_end    = if_else(N2S <= N2E, WSE, WSS),
    slope    = (W_end - W_start) / (N2_end - N2_start),
    # clip N2 coords to window
    x1_clip  = pmax(N2_start, win_start),
    x2_clip  = pmin(N2_end,   win_end),
    # project clipped N2 positions into W coords
    y1_clip  = W_start + (x1_clip - N2_start) * slope,
    y2_clip  = W_start + (x2_clip - N2_start) * slope
  ) %>%
  # drop any segments that vanish after clipping
  dplyr::filter(x1_clip < x2_clip)


conservedINV <- ggplot(threeINVnucmer_trimmed) +
  geom_segment(aes(x    = x1_clip / 1e6,xend = x2_clip / 1e6,y    = y1_clip / 1e6,yend = y2_clip / 1e6,color = contig), linewidth = 1) +
  geom_rect(data = threeINV %>% dplyr::slice_head(n = 1) %>% dplyr::select(chrom, pos, sv_length), aes(xmin = pos / 1e6, xmax = (pos + sv_length) / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(~strain, scales = "free") +
  coord_cartesian(xlim = c(win_start, win_end) / 1e6) +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)") +
  theme(
    legend.position  = "none",
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_blank(),
    panel.background = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_rect(fill = NA))
conservedINV
























# ======================================================================================================================================================================================== #
# Circos variation plot
# ======================================================================================================================================================================================== #
snps <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/misc/140WSs_biallelicSNPs.tsv", col_names = c("chrom","pos","ref","alt")) # don't currently have SNP calls for CGC1
merged_SV <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/jasmine_SVmerging/elegans/output/summarized_data.tsv")
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
geo_initial <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)
WSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/140_Ce_WSs.tsv", col_names = "strain") %>% dplyr::pull()
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") 
n2_genes_plt <- N2_gff %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::select(seqid,start,end,attributes) %>%
  dplyr::rename(chrom = seqid) %>% 
  dplyr::select(chrom, start, end) %>% 
  dplyr::filter(chrom != "MtDNA")

# Isolation site of each wild strain
geo <- geo_initial %>%
  dplyr::left_join(hawaii_islands, by = "isotype") %>%
  dplyr::mutate(geo = ifelse(geo == "Hawaii",collection_island_Hawaii,geo)) %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::filter(isotype %in% WSs)




# Merged SV plot
jasmine_plt <- ggplot(data = merged_SV %>% dplyr::filter(chrom != "MtDNA") %>% dplyr::select(chrom, pos, sv_type, number_svs_merged)) +
  geom_rect(aes(xmin = (pos / 1e6) - 0.001, xmax = (pos / 1e6) + 0.001, ymin = 0, ymax = number_svs_merged, fill = sv_type)) + 
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 14, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # strip.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,150)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "N2 genome position (Mb)", y = "Number of strains contributing to merged SV")
jasmine_plt

jasmine_plt_hist <- ggplot(data = merged_SV) +
  geom_histogram(aes(x = number_svs_merged, fill = sv_type), binwidth = 1, position = position_dodge(width = 0.8)) +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 14, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA)) +
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0,0)) +
  labs(x = "Number of strains contributing to merged SV", y = "Count")
jasmine_plt_hist

# ensuring merging of SV calls is correct
# merged_all <- merged_SV %>% dplyr::filter(number_svs_merged == "141") %>%
#   dplyr::mutate(sv_length = abs(sv_length)) %>%
#   dplyr::group_by(sv_type) %>%
#   dplyr::arrange(desc(sv_length)) %>%
#   dplyr::slice_head(n = 3) %>%
#   dplyr::ungroup()
# 
# 
nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain != "ECA396")
# 
# # INS found among all WSs
# ins_all <- ggplot(nucmer %>% dplyr::filter(N2_chr == "V")) +
#   geom_rect(data = merged_all %>% dplyr::filter(pos == "6175352"), aes(xmin = pos / 1e6, xmax = (pos + sv_length) / 1e6, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.5) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   facet_wrap(~strain) +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA)) +
#   coord_cartesian(xlim = c(6.17, 6.2)) +
#   labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
# ins_all
# 
# # DEL found among all WSs
# del_all <- ggplot(nucmer %>% dplyr::filter(N2_chr == "I")) +
#   geom_rect(data = merged_all %>% dplyr::filter(pos == "4760886"), aes(xmin = pos / 1e6, xmax = (pos + sv_length) / 1e6, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.5) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   facet_wrap(~strain) +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA)) +
#   coord_cartesian(xlim = c(4.760786, 4.762886)) +
#   labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
# del_all
# 
# 
# # INV found among 140 WSs - it is found in CGC1 - ECA1286 does not have the INV call
# inv_most <- merged_SV %>%
#   dplyr::filter(sv_type == "INV") %>% 
#   dplyr::mutate(sv_length = abs(sv_length)) %>%
#   dplyr::arrange(desc(number_svs_merged)) %>%
#   dplyr::slice_head(n = 3)
# 
# inv_140 <- ggplot(nucmer %>% dplyr::filter(N2_chr == "II", strain != "ECA1286")) +
#   geom_rect(data = inv_most %>% dplyr::filter(pos == "5408382"), aes(xmin = pos / 1e6, xmax = (pos + sv_length) / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   facet_wrap(~strain) +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA)) +
#   coord_cartesian(xlim = c(5.408182, 5.411382)) +
#   labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
# inv_140
# 
# eca1286 <- ggplot(nucmer %>% dplyr::filter(N2_chr == "II", strain == "ECA1286")) +
#   geom_rect(data = inv_most %>% dplyr::filter(pos == "5408382"), aes(xmin = pos / 1e6, xmax = (pos + sv_length) / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
#   geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
#   theme_bw() +
#   facet_wrap(~strain) +
#   theme(
#     legend.position = 'none',
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_text(size = 16, color = 'black', face = 'bold'),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(fill = NA)) +
#   coord_cartesian(xlim =c(5.408182, 5.411382)) +
#   labs(x = "N2 genome position (Mb)", y = "ECA1286 contig position (Mb)")
# eca1286
# 
# # need to visualize like how I visualized the 88 strain INV!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# strains <- nucmer %>% dplyr::select(strain) %>% dplyr::distinct() %>% dplyr::pull(strain) #%>% dplyr::filter(strain != "ECA1286") %>% dplyr::pull()
# 
# onefourty_INV <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
#   dplyr::filter(strain %in% strains) %>% dplyr::filter(N2_chr == "II") %>%
#   dplyr::filter(N2S < 5408182 & N2E < 5411382 & N2E > 5408182 | 
#                   N2S > 5408182 & N2S < 5411382 & N2E > 5411382 | 
#                   N2S > 5408182 & N2E < 5411382) #%>%
#   dplyr::group_by(strain, contig) %>%
#   dplyr::mutate(summedL2 = sum(L2)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(strain) %>%
#   dplyr::filter(summedL2 == max(summedL2)) %>%
#   dplyr::ungroup() %>% 
#   dplyr::mutate(slope = ifelse(N2S < 5408182 | N2E > 5411382, (WSS - WSE)/(N2S-N2E), NA)) %>%
#   dplyr::mutate(inv = ifelse(WSE < WSS, T, F)) %>%
#   dplyr::mutate(intercept = WSS - (slope * N2S)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
#   dplyr::mutate(new_start = ifelse(N2S < 5408182 & N2E < 5411382 & inv == F, ((slope * 5408182) + intercept), 
#                                    ifelse(N2S > 5408182 & N2E > 5411382  & inv == T, ((slope * 5408182) + intercept), WSS))) %>%
#   dplyr::mutate(new_end = ifelse(N2S > 5408182 & N2E > 5411382 & inv == F, ((slope * 5411382) + intercept), 
#                                  ifelse(N2S < 5408182 & N2E < 5411382 & inv == T, ((slope * 5408182) + intercept), WSE))) %>%
#   dplyr::select(-IDY,-L2,-L1,-LENR,-LENQ)
# 
# win_start <- 5408182
# win_end   <- 5411382
# 
# onefourty_INV_trimmed <- onefourty_INV %>%
#   # keep only alignments that overlap the window at all
#   dplyr::filter(N2E >= win_start, N2S <= win_end) %>%
#   # normalize so N2_start < N2_end, and match W coords to that direction
#   dplyr::mutate(
#     N2_start = pmin(N2S, N2E),
#     N2_end   = pmax(N2S, N2E),
#     W_start  = if_else(N2S <= N2E, WSS, WSE),
#     W_end    = if_else(N2S <= N2E, WSE, WSS),
#     slope    = (W_end - W_start) / (N2_end - N2_start),
#     # clip N2 coords to window
#     x1_clip  = pmax(N2_start, win_start),
#     x2_clip  = pmin(N2_end,   win_end),
#     # project clipped N2 positions into W coords
#     y1_clip  = W_start + (x1_clip - N2_start) * slope,
#     y2_clip  = W_start + (x2_clip - N2_start) * slope
#   ) %>%
#   # drop any segments that vanish after clipping
#   dplyr::filter(x1_clip < x2_clip)
# 
# 
# conservedINV <- ggplot(onefourty_INV_trimmed) +
#   geom_segment(aes(x    = x1_clip / 1e6,xend = x2_clip / 1e6,y    = y1_clip / 1e6,yend = y2_clip / 1e6,color = contig), linewidth = 1) +
#   geom_rect(data = inv_most %>% dplyr::filter(pos == "5408382"), aes(xmin = pos / 1e6, xmax = (pos + sv_length) / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
#   facet_wrap(~strain, scales = "free") +
#   coord_cartesian(xlim = c(win_start, win_end) / 1e6) +
#   labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)") +
#   theme(
#     legend.position  = "none",
#     axis.text        = element_blank(),
#     axis.ticks       = element_blank(),
#     axis.title       = element_blank(),
#     panel.background = element_blank(),
#     panel.grid       = element_blank(),
#     panel.border     = element_rect(fill = NA))
# conservedINV





# N2 genes plot
genes_plt <- ggplot(n2_genes_plt %>% dplyr::filter(chrom != "MtDNA")) + 
  geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymin = 0, ymax = 1), fill = "black") +
  facet_wrap(~chrom, scales = "free_x", nrow = 1) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 
genes_plt


# Calculating snp count per kb
snps <- snps %>%
  dplyr::mutate(bin = (pos %/% 1000) * 1000) %>%  # floor to nearest 1000
  dplyr::count(chrom, bin, name = "variant_count") %>%
  dplyr::arrange(chrom, bin)

# SNP density plot
snps_plt <- ggplot(snps) + 
  geom_point(aes(x = bin, y = variant_count), color = '#DB6333', alpha = 0.7) +
  facet_wrap( ~chrom, nrow = 1, scales = "free_x") + 
  # geom_smooth(aes(x = bin, y = variant_count), method = "loess", se = TRUE, color = "lightblue") +
  # ylab("Variants per kb") + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12, color = 'black'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) 
snps_plt 



# Calculating HDR frequency
bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

strains <- nucmer %>% dplyr::select(strain) %>% dplyr::distinct() %>% dplyr::filter(strain != "CGC1") %>% dplyr::pull()
hdrs <- hdrs %>% dplyr::filter(strain %in% strains)
hdrs_dt <- as.data.table(hdrs)

setkey(bins_dt, chrom, start, end)
setkey(hdrs_dt, chrom, start, end)

overlaps <- data.table::foverlaps(hdrs_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/140) #change me to number of isotypes

hdr_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

hdr_freq <- ggplot(hdr_bin_plt) +
  geom_point(aes(x = middle, y = freq), color = 'gray') +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12, color = 'black'),
    panel.grid = element_blank(),
    # panel.grid.major= element_line(color = 'gray80'),
    # panel.grid.major.x = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) +
  coord_cartesian(ylim = c(0,1))
hdr_freq





# COLLAPSING HDRS AMONG 140 WSs
all_regions <- hdrs %>%
  dplyr::rename(START = start, CHROM = chrom, END = end) %>%
  dplyr::arrange(CHROM,START) %>%
  dplyr::group_split(CHROM)

strain_count <- hdrs %>% dplyr::distinct(strain, .keep_all = T)
print(nrow(strain_count)) # SHOULD BE 140

# Collapsing all HDRs
getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>%
        dplyr::arrange(CHROM,START) %>%
        dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))

      #print(nrow(checkIntersect %>% dplyr::filter(check==T)))

      if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
        print("NO MORE INTERSECTS")
        k=0
      } else {

        temp <- checkIntersect %>%
          dplyr::mutate(gid=data.table::rleid(check)) %>%
          dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))

        collapse <- temp %>%
          dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
          dplyr::group_by(gid) %>%
          dplyr::mutate(newStart=min(START)) %>%
          dplyr::mutate(newEnd=max(END)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(START=newStart,END=newEnd) %>%
          dplyr::select(-newEnd,-newStart)

        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)

        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)

        j=j+1
      }
    }
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

HDR_collapse_master <- getRegFreq(all_regions)

all_collapsed <- plyr::ldply(HDR_collapse_master, data.frame) %>%
  dplyr::select(-strain)

colnames(all_collapsed) <- c("chrom","start","end")



# Calculating DEL frequency
bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

del_calls <- filt_calls %>% dplyr::filter(sv_type == "DEL", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
del_calls_dt <- as.data.table(del_calls)

setkey(bins_dt, chrom, start, end)
setkey(del_calls_dt, chrom, start, end)

overlaps <- data.table::foverlaps(del_calls_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/141) #change me to number of isotypes

del_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

del_freq <- ggplot(del_bin_plt) +
  geom_rect(data = all_collapsed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'gray', alpha = 0.7) +
  geom_point(aes(x = middle, y = freq), color = 'red') +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12, color = 'black'),
    # panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0))
del_freq


# Calculating INS frequency
bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

ins_calls <- filt_calls %>% dplyr::filter(sv_type == "INS", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
ins_calls_dt <- as.data.table(ins_calls)

setkey(bins_dt, chrom, start, end)
setkey(ins_calls_dt, chrom, start, end)

overlaps <- data.table::foverlaps(ins_calls_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/141) #change me to number of isotypes

ins_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

ins_freq <- ggplot(ins_bin_plt) +
  geom_rect(data = all_collapsed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'gray', alpha = 0.7) +
  geom_point(aes(x = middle, y = freq), color = 'blue') +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12, color = 'black'),
    # panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0))
ins_freq



# Calculating INV frequency
bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
bins_dt <- as.data.table(bins)
bins_dt[, id := .I]  # optional: keep track of bins

inv_calls <- filt_calls %>% dplyr::filter(sv_type == "INV", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)
inv_calls_dt <- as.data.table(inv_calls)

setkey(bins_dt, chrom, start, end)
setkey(inv_calls_dt, chrom, start, end)

overlaps <- data.table::foverlaps(inv_calls_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("chrom", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq = n_strains/141) #change me to number of isotypes

inv_bin_plt <- bins_wFreq %>% dplyr::mutate(middle = (end + start) / 2)

inv_freq <- ggplot(inv_bin_plt) +
  geom_rect(data = all_collapsed, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), fill = 'gray', alpha = 0.7) +
  geom_point(aes(x = middle, y = freq), color = 'gold1') +
  facet_wrap(~chrom, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 12, color = 'black'),
    # panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) +
  coord_cartesian(ylim = c(0,1.01)) +
  scale_y_continuous(expand = c(0,0))
inv_freq


all_three <- cowplot::plot_grid(
  del_freq, ins_freq, inv_freq,
  nrow = 3,
  align = "v",
  rel_heights = c(1,1,1))
all_three

# =============================================================== #
# Enrichment of SVs in HDRs #
# =============================================================== #
library(GenomicRanges)
library(IRanges)

chrom_sizes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/N2.WS283.cleaned.fa.fai", col_names = c("chrom","start","end"))

# DEL enrichment in HDRs
hdr_bins <- all_collapsed
del_calls <- filt_calls %>% dplyr::filter(sv_type == "DEL", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)

# helper: tibble (BED-like, 0-based) -> GRanges (1-based)
to_gr <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start + 1, end = df$end)
  )
}

genome_gr <- to_gr(chrom_sizes)          # whole genome space
hdr_gr    <- to_gr(hdr_bins)  
del_gr    <- to_gr(del_calls)            # DEL intervals

# de-duplicate identical DELs (non-rare HDRs) because we are assessing genomic loci affected by HDRs in or outside of HDRs, not the frequency that a genomic
# locus is affected by a DEL
# Keeping duplicated DEL calls - "Among all DEL calls across 140 strain, are calls more frequent per bp in HDRs than outside?"
# De-duplicating DEL calls - "Are genomic loci that harbor deletions enriched in HDRs vs outside?"
del_gr_uniq <- unique(del_gr)

# count DEL events overlapping HDR
in_hdr <- countOverlaps(del_gr_uniq, hdr_gr) > 0
n_hdr <- sum(in_hdr) # number of DELs that overlap with HDRs
n_out <- length(del_gr_uniq) - n_hdr # number of DELs that do NOT overlap with HDRs

# bp denominators
hdr_bp    <- sum(width(hdr_gr)) # summed genome size of HDRs
genome_bp <- sum(width(genome_gr)) # genome size
nonhdr_bp <- genome_bp - hdr_bp # summed genome size outside of HDRs

# densities and fold enrichment
density_hdr <- n_hdr / hdr_bp # normalizing by bps - # of DELs in HDRs / # bps in HDRs
density_out <- n_out / nonhdr_bp

fold_enrichment <- density_hdr / density_out

# report
list(
  total_DEL = length(del_gr_uniq),
  DEL_in_HDR = n_hdr,
  DEL_out_HDR = n_out,
  HDR_bp = hdr_bp,
  nonHDR_bp = nonhdr_bp,
  density_HDR = density_hdr,
  density_nonHDR = density_out,
  fold_enrichment = fold_enrichment)


# INS enrichment in HDRs
hdr_bins <- all_collapsed
ins_calls <- filt_calls %>% dplyr::filter(sv_type == "INS", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)

to_gr <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start + 1, end = df$end)
  )
}

genome_gr <- to_gr(chrom_sizes)          # whole genome space
hdr_gr    <- to_gr(hdr_bins)  
ins_gr    <- to_gr(ins_calls)            # ins intervals


ins_gr_uniq <- unique(ins_gr)

in_hdr <- countOverlaps(ins_gr_uniq, hdr_gr) > 0
n_hdr <- sum(in_hdr) # number of INSs that overlap with HDRs
n_out <- length(ins_gr_uniq) - n_hdr # number of INSs that do NOT overlap with HDRs

hdr_bp    <- sum(width(hdr_gr)) # summed genome size of HDRs
genome_bp <- sum(width(genome_gr)) # genome size
nonhdr_bp <- genome_bp - hdr_bp # summed genome size outside of HDRs

density_hdr <- n_hdr / hdr_bp # normalizing by bps - # of INSs in HDRs / # bps in HDRs
density_out <- n_out / nonhdr_bp

fold_enrichment <- density_hdr / density_out

# report
list(
  total_ins = length(ins_gr_uniq),
  ins_in_HDR = n_hdr,
  ins_out_HDR = n_out,
  HDR_bp = hdr_bp,
  nonHDR_bp = nonhdr_bp,
  density_HDR = density_hdr,
  density_nonHDR = density_out,
  fold_enrichment = fold_enrichment)


# INV enrichment in HDRs
hdr_bins <- all_collapsed
INV_calls <- filt_calls %>% dplyr::filter(sv_type == "INV", strain != "CGC1") %>% dplyr::mutate(end = pos + sv_length) %>% dplyr::rename(start = pos) %>% dplyr::select(chrom,start,end,strain)

to_gr <- function(df) {
  GRanges(
    seqnames = df$chrom,
    ranges   = IRanges(start = df$start + 1, end = df$end)
  )
}

genome_gr <- to_gr(chrom_sizes)          # whole genome space
hdr_gr    <- to_gr(hdr_bins)  
INV_gr    <- to_gr(INV_calls)            # ins intervals


INV_gr_uniq <- unique(INV_gr)

in_hdr <- countOverlaps(INV_gr_uniq, hdr_gr) > 0
n_hdr <- sum(in_hdr) # number of INSs that overlap with HDRs
n_out <- length(INV_gr_uniq) - n_hdr # number of INSs that do NOT overlap with HDRs

hdr_bp    <- sum(width(hdr_gr)) # summed genome size of HDRs
genome_bp <- sum(width(genome_gr)) # genome size
nonhdr_bp <- genome_bp - hdr_bp # summed genome size outside of HDRs

density_hdr <- n_hdr / hdr_bp # normalizing by bps - # of INSs in HDRs / # bps in HDRs
density_out <- n_out / nonhdr_bp

fold_enrichment <- density_hdr / density_out

# report
list(
  total_ins = length(INV_gr_uniq),
  INV_in_HDR = n_hdr,
  INV_out_HDR = n_out,
  HDR_bp = hdr_bp,
  nonHDR_bp = nonhdr_bp,
  density_HDR = density_hdr,
  density_nonHDR = density_out,
  fold_enrichment = fold_enrichment)




# ###########################
# Circos plot - example
# ###########################
# col_fun = colorRamp2(c(-2, 0, 2), c("pink", "purple", "red"))
# circlize_plot = function() {
#   set.seed(12345)
#   sectors = letters[1:12]
#   circos.initialize(sectors, xlim = c(0, 1))
#   circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#     circos.points(runif(20), runif(20), cex = 0.5, pch = 16, col = 2)
#     circos.points(runif(20), runif(20), cex = 0.5, pch = 16, col = 3)
#   })
#   circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
#     circos.lines(sort(runif(20)), runif(20), col = 4)
#     circos.lines(sort(runif(20)), runif(20), col = 5)
#   })
#   
#   for(i in 1:12) {
#     circos.link(sample(sectors, 1), sort(runif(10))[1:2], 
#                 sample(sectors, 1), sort(runif(10))[1:2],
#                 col = add_transparency(col_fun(rnorm(1))))
#   }
#   circos.clear()
# }
# circlize_plot()





## =========================
## CIRCLIZE PLOT (ALL-IN-ONE)
## Outer chromosome ring with genes plotted + inner tracks:
##   1) SNP density (points)
##   2) HDR frequency (points)
##   3) DEL frequency (points)
##   4) INS frequency (points)
##   5) INV frequency (points)
## Includes legend (ComplexHeatmap::Legend)
## =========================

bin_size <- 1000
chr_order <- c("I","II","III","IV","V","X")

## Path to genome .fai (recommended; provides chromosome lengths)
## If you don't have this, see fallback block below.
fai_path <- "genome.fa.fai"

## ---- INPUT OBJECTS ASSUMED TO EXIST (from your pipeline) ----
## n2_genes_plt: data.frame with chrom, start, end (bp)
## snps: tibble/data.frame with chrom, pos, ref, alt (pos in bp)
## hdrs: tibble/data.frame with chrom, start, end, strain (BED-like)
## filt_calls: tibble/data.frame with chrom, pos, sv_length, sv_type, strain
## nucmer: tibble/data.frame with strain (used to define subset of strains)
##
## If any are not loaded, load them before running this block.

## =========================
## 1) CHROM SIZES
## =========================
chrom_sizes <- readr::read_tsv(fai_path, col_names = c("chrom","len","x1","x2","x3"), show_col_types = FALSE) %>%
  dplyr::select(chrom, len) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(chrom = factor(chrom, levels = chr_order)) %>%
  dplyr::arrange(chrom)

## Fallback (if you don't have .fai):
## Uncomment to infer lengths from bins/genes if needed.
# chrom_sizes <- bind_rows(
#   snps %>% group_by(chrom) %>% summarise(len = max(pos, na.rm=TRUE), .groups="drop"),
#   n2_genes_plt %>% group_by(chrom) %>% summarise(len = max(end, na.rm=TRUE), .groups="drop")
# ) %>% group_by(chrom) %>% summarise(len = max(len), .groups="drop") %>%
#   filter(chrom != "MtDNA") %>%
#   mutate(chrom = factor(chrom, levels = chr_order)) %>% arrange(chrom)


## =========================
## 3) SNP DENSITY PER BIN (counts per bin_size)
## =========================
snps_binned <- snps %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(bin = (pos %/% bin_size) * bin_size) %>%
  dplyr::count(chrom, bin, name = "variant_count") %>%
  dplyr::mutate(pos = bin + bin_size/2) %>%
  dplyr::select(chrom, pos, value = variant_count)

## =========================
## 4) HDR FREQUENCY PER BIN
##    (subset strains from nucmer, excluding CGC1, as you did)
## =========================
strains_hdr <- nucmer %>%
  dplyr::distinct(strain) %>%
  dplyr::filter(strain != "CGC1") %>%
  dplyr::pull(strain)
n_hdr <- length(strains_hdr)

hdrs_sub <- hdrs %>% dplyr::filter(strain %in% strains_hdr, chrom != "MtDNA")
hdrs_dt  <- as.data.table(hdrs_sub)
setkey(hdrs_dt, chrom, start, end)

ov_hdr <- foverlaps(hdrs_dt, bins_dt, nomatch = 0L)

hdr_counts <- ov_hdr[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end, id)]
hdr_w <- merge(bins_dt, hdr_counts, by = c("chrom","start","end","id"), all.x = TRUE)
hdr_w[is.na(n_strains), n_strains := 0L]
hdr_w[, value := n_strains / n_hdr]
hdr_track <- as.data.frame(hdr_w) %>%
  dplyr::mutate(pos = (start + end)/2) %>%
  dplyr::select(chrom, pos, value)

## =========================
## 5) SV FREQUENCY PER BIN (DEL/INS/INV)
## =========================
sv_freq_track <- function(filt_calls, svtype, bins_dt, denom_n) {
  calls <- filt_calls %>%
    dplyr::filter(sv_type == svtype, chrom != "MtDNA") %>%
    dplyr::mutate(start = pos,
                  end   = pos + sv_length) %>%
    dplyr::select(chrom, start, end, strain)
  
  calls_dt <- as.data.table(calls)
  setkey(calls_dt, chrom, start, end)
  
  ov <- foverlaps(calls_dt, bins_dt, nomatch = 0L)
  counts <- ov[, .(n_strains = uniqueN(strain)), by = .(chrom, start, end, id)]
  
  w <- merge(bins_dt, counts, by = c("chrom","start","end","id"), all.x = TRUE)
  w[is.na(n_strains), n_strains := 0L]
  w[, value := n_strains / denom_n]
  
  as.data.frame(w) %>%
    dplyr::mutate(pos = (start + end)/2) %>%
    dplyr::select(chrom, pos, value)
}

## You used /141 for SV frequencies; keep your convention:
n_sv <- 141

del_track <- sv_freq_track(filt_calls, "DEL", bins_dt, n_sv)
ins_track <- sv_freq_track(filt_calls, "INS", bins_dt, n_sv)
inv_track <- sv_freq_track(filt_calls, "INV", bins_dt, n_sv)

## =========================
## 6) CIRCLIZE PLOT
## =========================
add_point_track <- function(df, col, ylim, track_height = 0.10, label = NULL) {
  circos.trackPlotRegion(
    ylim = ylim,
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- df[df$chrom == chr, , drop = FALSE]
      if (nrow(d) == 0) return()
      circos.points(d$pos, d$value, pch = 16, cex = 0.25, col = col)
    }
  )
  if (!is.null(label)) {
    # Put track label near the first chromosome sector
    circos.text(
      x = 0, y = mean(ylim), labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}

## Open a device if desired (e.g., PDF):
## pdf("circos_tracks.pdf", width=8, height=8)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = c(rep(2, length(chr_order)-1), 8),
  track.margin = c(0.002, 0.002),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(
  factors = chrom_sizes$chrom,
  xlim = cbind(rep(0, nrow(chrom_sizes)), chrom_sizes$len)
)

## Outer chromosome ring (ideogram-like) - I NEED TO ADD GENE MODELS
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.065,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    
    circos.rect(xlim[1], 0, xlim[2], 1, col = "grey90", border = "white")
    circos.text(CELL_META$xcenter, 0.55, chr, facing = "bending.inside", cex = 0.9, font = 2)
    
    ## Axis ticks (Mb)
    circos.axis(
      h = "top",
      major.at = seq(0, xlim[2], by = 5e6),
      labels = seq(0, xlim[2], by = 5e6) / 1e6,
      labels.cex = 0.5,
      major.tick.length = 0.02,
      minor.ticks = 4
    )
  }
)

## Tracks
snp_ylim <- c(0, max(snps_binned$value, na.rm = TRUE))
add_point_track(snps_binned, col = "#DB6333", ylim = snp_ylim, track_height = 0.12, label = "SNPs/kb")
add_point_track(hdr_track,  col = "grey40",  ylim = c(0,1),   track_height = 0.10, label = "HDR freq")
add_point_track(del_track,  col = "red",     ylim = c(0,1),   track_height = 0.08, label = "DEL freq")
add_point_track(ins_track,  col = "blue",    ylim = c(0,1),   track_height = 0.08, label = "INS freq")
add_point_track(inv_track,  col = "gold",    ylim = c(0,1),   track_height = 0.08, label = "INV freq")

## Legend (place top-right)
lgd <- Legend(
  labels = c("SNP density", "HDR freq", "DEL freq", "INS freq", "INV freq"),
  type = "points",
  pch = 16,
  legend_gp = gpar(col = c("#DB6333", "grey40", "red", "blue", "gold"))
)
draw(lgd, x = unit(0.88, "npc"), y = unit(0.88, "npc"))

circos.clear()

































# ============================================== # 
# SVs in coding regions
# ============================================== # 
MAF_thresh <- round(0.05 * 141)

maf_filt <- merged_SV %>% 
  dplyr::select(chrom,pos,sv_length,sv_type,number_svs_merged) %>%
  dplyr::filter(number_svs_merged > MAF_thresh) %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::mutate(end = pos + sv_length) %>% 
  dplyr::rename(start = pos) %>%
  dplyr::select(chrom, start, end, sv_type) %>%
  dplyr::mutate(overlap = F)

svs_dt <- as.data.table(maf_filt)
n2_genes_dt <- as.data.table(n2_genes_plt)

setkey(svs_dt, chrom, start, end)
setkey(n2_genes_dt, chrom, start, end)

svs_inCodingRegions <- data.table::foverlaps(x = svs_dt, y = n2_genes_dt, type = "any") %>% dplyr::filter(!is.na(start)) %>% dplyr::mutate(overlap = T)


# Ensuring that the foverlaps command worked correctly
test <- svs_inCodingRegions %>% dplyr::filter(start > 1600000 & end < 1700000) %>% dplyr::select(chrom,start,end) %>% dplyr::distinct(chrom,start,end)

check <- ggplot(svs_inCodingRegions %>% dplyr::filter(start > 1600000 & end < 1700000)) +
  geom_rect(data = test, aes(xmin = start / 1e6, xmax = end /1e6, ymin = 0, ymax = 0.99, fill = "N2_genes")) +
  geom_rect(aes(xmin = i.start / 1e6, xmax = i.end /1e6, ymin = 1.01, ymax = 2, fill = "SVs")) +
  scale_fill_manual(values = c("N2_genes" = "forestgreen", "SVs" = "purple")) + 
  facet_wrap(~chrom, nrow = 2, scales = "free_x") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, color = 'black'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black")) 
check


overlap <- svs_inCodingRegions %>% dplyr::select(chrom, i.start,i.end, sv_type, overlap) %>% dplyr::rename(start = i.start, end = i.end) %>% dplyr::distinct(chrom,start,end,sv_type, .keep_all = T)

final_stats <- maf_filt %>% 
  dplyr::left_join(overlap, by = c("chrom", "start", "end", "sv_type")) %>% 
  dplyr::mutate(region = ifelse(is.na(overlap.y),'non-PC_region','overlaps_PCgene')) %>%
  dplyr::select(chrom,start,end,sv_type,region) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::mutate(total_sv_type = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sv_type,region) %>%
  dplyr::mutate(region_count = n()) %>%
  dplyr::ungroup()

plt_stats <- final_stats %>% dplyr::select(sv_type,region,total_sv_type,region_count) %>%
  dplyr::mutate(proportion = (region_count / total_sv_type) * 100) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS","DEL","INV")))


ggplot() +
  geom_bar(data = plt_stats, aes(x = sv_type, y = proportion, fill = region), stat = "identity") +
  geom_text(data = plt_stats, aes(x = sv_type, y = proportion, label = region_count, group = region), position = position_stack(vjust = 0.5),color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("overlaps_PCgene" = "red", "non-PC_region" = "gray40")) +
  labs(y = "Proportion (%)", fill = "Region") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 12, face = 'bold'),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))


overlapped_n2_genes <- svs_inCodingRegions %>% dplyr::select(chrom, start, end) %>% dplyr::distinct() %>% dplyr::mutate(class = "overlapped")
prop_n2_genes <- n2_genes_plt %>% dplyr::left_join(overlapped_n2_genes, by = c("chrom","start","end")) %>% dplyr::mutate(class = ifelse(is.na(class),"no_SV_overlap",class))


# .... load in variables "count", "ortho_genes_dd", and "private_OGs" from orthogroup_vis.R ............................... #
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))

strainCol <- colnames(ortho_genes_dd)
ugh <- gsub(".20251012.inbred.blobFiltered.softMasked.braker.longestIso.protein","", strainCol)
ugh2 <- gsub(".20251014.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh)
ugh3 <- gsub(".20251124.inbred.blobFiltered.softMasked.braker.longestIso.protein","", ugh2)
ugh4 <- gsub(".20251012.inbred.onlyONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh3)
ugh5 <- gsub(".Nov2025.softMasked.braker.longest.protein","", ugh4)
ugh6 <- gsub(".20251012.inbred.withONT.blobFiltered.softMasked.braker.longestIso.protein","", ugh5)
strainCol_c2 <- gsub("c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein","N2", ugh6)
colnames(ortho_genes_dd) <- strainCol_c2

ortho_count <- ortho_genes_dd

strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations_pre <- ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))


private_OGs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Orthogroups/Orthogroups_UnassignedGenes.tsv") %>%
  dplyr::filter(!grepl("MTCE",c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein))

colnames(private_OGs) <- strainCol_c2

private_cols <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup")]

private_ortho_count <- private_OGs
for (i in 1:length(private_cols)) {
  print(paste0(i, " out of ", length(private_cols)))
  temp_colname <- paste0(private_cols[i], "_count")
  
  private_ortho_count <- private_ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := ifelse(is.na(!!sym(private_cols[i])), NA, 1))
}

all_relations_private <- private_ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))

all_relations <- all_relations_pre %>%
  dplyr::bind_rows(all_relations_private)

private_freq = (1/(length(strainCol_c2_u)))

count <- all_relations %>%
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
  )
# ............................................................................................................................ #

n2_gene <- N2_gff %>%
  dplyr::filter(type == "gene") %>%
  tidyr::separate(attributes, into = c("blah","first"), sep = ";") %>%
  tidyr::separate(first, into = c('blah2','tran'), sep = ',') %>%
  dplyr::mutate(blah2 = gsub("Alias=","",blah2)) %>%
  dplyr::mutate(final_tran = ifelse(is.na(tran), blah2, tran)) %>%
  dplyr::select(-blah,-blah2,-tran) %>%
  dplyr::rename(tran = final_tran) %>%
  dplyr::mutate(tran = paste0("transcript_",tran))

n2_table <- ortho_genes_dd %>%
  dplyr::bind_rows(private_OGs) %>% 
  dplyr::select(Orthogroup,N2)

ortho_count_wCoord <- count %>%
  dplyr::left_join(n2_table, by = "Orthogroup") %>%
  dplyr::select(freq, class, N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ",\\s*") %>% # splitting rows so each gene is on a row and it retains is gene set classification
  ### REMOVE .* for isoform in order to left join!
  dplyr::mutate(N2 = sub("\\.[^.]*$", "", N2)) %>% # removing the last period in the transcript name and the isoform number
  dplyr::mutate(N2 = sub("[A-Za-z]$", "", N2)) %>% # removing any letters that are in the last position in a string
  dplyr::left_join(n2_gene, by = c("N2" = "tran")) %>%
  dplyr::select(freq,class,N2,seqid,start,end) %>%
  dplyr::distinct(N2, seqid, start, end, .keep_all = T) %>%
  dplyr::group_by(N2) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

plt_data <- ortho_count_wCoord %>%
  dplyr::filter(seqid != "MtDNA") %>%
  dplyr::mutate(mid_mb = (start + end) / 2 / 1e6)

n2_genes_geneSet <- plt_data %>% dplyr::select(seqid,start,end,class) %>% dplyr::rename(chrom = seqid, gene_set = class) %>%
  dplyr::left_join(prop_n2_genes, by = c("chrom","start","end")) %>% 
  dplyr::mutate(final_class = ifelse(class == "overlapped",gene_set,class)) %>% 
  dplyr::select(final_class) %>% 
  dplyr::group_by(final_class) %>% 
  dplyr::mutate(total_class = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total = dplyr::n()) %>%
  dplyr::distinct() %>%
  dplyr::mutate(prop = total_class / total * 100) %>% 
  dplyr::mutate(final_class = ifelse(final_class == "no_SV_overlap","none",final_class)) %>%
  dplyr::mutate(final_class = factor(final_class, levels = c("none","private","accessory","core")))

ggplot() +
  geom_bar(data = n2_genes_geneSet, aes(x = "N2 genes", y = prop, fill = final_class), stat = "identity") +
  geom_text(data = n2_genes_geneSet, aes(x = "N2 genes", y = prop, label = total_class, group = final_class), position = position_stack(vjust = 0.5),color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("core" = "green4", "none" = "gray40", "accessory" = "#DB6333", "private" = "magenta3")) +
  labs(y = "Proportion (%)", fill = "SV overlap with N2 genes") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 12, face = 'bold'),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))








### Overlap with an N2 gene or 2kb upstream (encompases promoter region)
twokb_n2_genes_plt <- n2_genes_plt %>% dplyr::mutate(start = start - 2000)
svs_dt <- as.data.table(maf_filt)
n2_genes_dt <- as.data.table(twokb_n2_genes_plt)

setkey(svs_dt, chrom, start, end)
setkey(n2_genes_dt, chrom, start, end)

svs_inCodingRegions <- data.table::foverlaps(x = svs_dt, y = n2_genes_dt, type = "any") %>% dplyr::filter(!is.na(start)) %>% dplyr::mutate(overlap = T)

overlap <- svs_inCodingRegions %>% dplyr::select(chrom, i.start,i.end, sv_type, overlap) %>% dplyr::rename(start = i.start, end = i.end) %>% dplyr::distinct(chrom,start,end,sv_type, .keep_all = T)

final_stats <- maf_filt %>% 
  dplyr::left_join(overlap, by = c("chrom", "start", "end", "sv_type")) %>% 
  dplyr::mutate(region = ifelse(is.na(overlap.y),'non-PC_region','overlaps_PCgene_(plus2kbupstream)')) %>%
  dplyr::select(chrom,start,end,sv_type,region) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::mutate(total_sv_type = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(sv_type,region) %>%
  dplyr::mutate(region_count = n()) %>%
  dplyr::ungroup()

plt_stats <- final_stats %>% dplyr::select(sv_type,region,total_sv_type,region_count) %>%
  dplyr::mutate(proportion = (region_count / total_sv_type) * 100) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS","DEL","INV")))


ggplot() +
  geom_bar(data = plt_stats, aes(x = sv_type, y = proportion, fill = region), stat = "identity") +
  geom_text(data = plt_stats, aes(x = sv_type, y = proportion, label = region_count, group = region), position = position_stack(vjust = 0.5), color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("overlaps_PCgene_(plus2kbupstream)" = "firebrick", "non-PC_region" = "gray60")) +
  labs(y = "Proportion (%)", fill = "Region") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 12, face = 'bold'),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))


# ============================================== # 
# SVs overlapping with specific N2 gene classes
# ============================================== # 
common_SVs <- maf_filt %>% 
  dplyr::select(-overlap)

common_SVs <- merged_SV %>% 
  dplyr::select(chrom,pos,sv_length,sv_type,number_svs_merged) %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::mutate(end = pos + sv_length) %>% 
  dplyr::rename(start = pos) %>%
  dplyr::select(chrom, start, end, sv_type) 


N2_tranGene <- N2_gff %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::mutate(attributes = gsub("ID=transcript:","",attributes), attributes = gsub("Parent=gene:","",attributes)) %>%
  tidyr::separate_wider_delim(attributes, delim = ";",names = c("tran", "N2", "rest"), too_many = "merge") %>%
  dplyr::select(tran, N2, -rest) %>%
  dplyr::mutate(tran = paste0("transcript:",tran))

cleaned_N2_gff <- N2_gff %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::select(seqid,start,end,attributes) %>%
  dplyr::rename(chrom = seqid) %>% 
  dplyr::filter(chrom != "MtDNA") %>% 
  dplyr::rename(N2 = attributes)

ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/ipr/output/N2_IPR_allApps_20251019.tsv", 
                       col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::left_join(N2_tranGene, by = 'tran') %>%
  dplyr::select(-tran) %>%
  dplyr::select(N2,signature_accession,signature_description,IPR_accession,IPR_description,GO)

ipr_gpcrs <- ipr %>% dplyr::filter(grepl("7TM", IPR_description)) %>% dplyr::select(N2, IPR_description) %>% dplyr::distinct(N2, .keep_all = T)
gpcrs <- cleaned_N2_gff %>%
  dplyr::left_join(ipr_gpcrs, by = "N2") %>% 
  dplyr::select(chrom, start, end, N2, IPR_description) %>% 
  dplyr::filter(!is.na(IPR_description))

gpcrs_list <- gpcrs %>% dplyr::pull(N2)
total_gpcr <- length(gpcrs_list)

# gpc_test <- ipr %>% dplyr::filter(grepl("7TM", IPR_description)) %>%
#   dplyr::group_by(N2) %>%
#   dplyr::mutate(annotation_count = n()) %>%
#   dplyr::ungroup()
#  
# multiple <- gpc_test %>% dplyr::filter(annotation_count > 1) %>% dplyr::distinct(N2)
#  
# at_least_one_anno <- gpc_test %>% dplyr::distinct(N2)

svs_dt <- as.data.table(common_SVs)
gpcrs_dt <- as.data.table(gpcrs)

setkey(svs_dt, chrom, start, end)
setkey(gpcrs_dt, chrom, start, end)

svs_inCodingRegions <- data.table::foverlaps(x = svs_dt, y = gpcrs_dt, type = "any") %>% dplyr::filter(!is.na(start)) %>% dplyr::mutate(overlap = T)

overlap <- svs_inCodingRegions %>% dplyr::select(chrom, start, end, N2, sv_type, overlap) %>% dplyr::distinct(chrom,start,end,N2,sv_type, .keep_all = T)

plt_stats_juvenile <- overlap %>%
  dplyr::filter(N2 %in% gpcrs_list) %>%
  dplyr::distinct(N2, sv_type) %>%          
  dplyr::count(sv_type, name = "n_genes") %>%
  dplyr::mutate(total_gpcr = total_gpcr, proportion = 100 * n_genes / total_gpcr) %>%
  dplyr::mutate(region = T)

plt_stats <- plt_stats_juvenile %>%
  tidyr::complete(sv_type, fill = list(n_genes = 0, region = "TRUE", total_gpcr = total_gpcr)) %>%
  dplyr::mutate(n_false = total_gpcr - n_genes) %>%
  select(sv_type, n_true = n_genes, n_false, total_gpcr) %>%
  pivot_longer(c(n_true, n_false), names_to = "region", values_to = "n_genes") %>%
  dplyr::mutate(region = recode(region, n_true = "TRUE", n_false = "FALSE"),proportion = 100 * n_genes / total_gpcr)

ggplot() +
  geom_bar(data = plt_stats, aes(x = sv_type, y = proportion, fill = region), stat = "identity") +
  geom_text(data = plt_stats, aes(x = sv_type, y = proportion, label = n_genes, group = region), position = position_stack(vjust = 0.5), color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "olivedrab", "FALSE" = "gray60")) +
  labs(y = "Proportion of N2 genes (%)", fill = "SV overlap with GPCR") +
  theme(panel.border = element_rect(color = 'black', fill = NA),
        panel.background = element_blank(),
        panel.grid.major= element_line(color = 'gray80'),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(color = 'black', size = 12, face = 'bold'),
        axis.title.y = element_text(size = 14, color = 'black', face = 'bold')) +
  guides(color = "none") + # to get rid of legend for the horizontal lines
  scale_y_continuous(expand = c(0,0))

