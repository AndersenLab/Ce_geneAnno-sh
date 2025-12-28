library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
install.packages("")

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
  coord_cartesian(xlims = c(win_start, win_end) / 1e6) +
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










































# chimeric??? 
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

























# ======================================================================================================================================================================================== #
# Circos variation plot
# ======================================================================================================================================================================================== #
snps <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/misc/140WSs_biallelicSNPs.tsv", col_names = c("chrom","pos","ref","alt"))
merged_SV <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/jasmine_SVmerging/output/summarized_data.tsv")
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
geo_initial <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)
WSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/140_Ce_WSs.tsv", col_names = "strain") %>% dplyr::pull()
n2_genes_plt <- n2_genes %>%
  dplyr::rename(chrom = seqid) %>% 
  dplyr::select(chrom, start, end)

# Isolation site of each wild strain
geo <- geo_initial %>%
  dplyr::left_join(hawaii_islands, by = "isotype") %>%
  dplyr::mutate(geo = ifelse(geo == "Hawaii",collection_island_Hawaii,geo)) %>%
  dplyr::select(isotype, lat, long, geo) %>%
  dplyr::filter(isotype %in% WSs)


# Calculating snp count per kb
snps <- snps %>%
  dplyr::mutate(bin = (pos %/% 1000) * 1000) %>%  # floor to nearest 1000
  dplyr::count(chrom, bin, name = "variant_count") %>%
  dplyr::arrange(chrom, bin)


# Calculating HDR frequency
bins <- snps %>% dplyr::select(chrom,bin) %>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::
bins_dt <- as.data.table(bins)
setnames(bins_dt, c("binStart", "binEnd"), c("start", "end"))
bins_dt[, id := .I]  # optional: keep track of bins

hdrs_dt <- as.data.table(hdrs)
setnames(hdrs_dt, c("minStart", "maxEnd"), c("start", "end"))

setkey(bins_dt, CHROM, start, end)
setkey(hdrs_dt, CHROM, start, end)

overlaps <- foverlaps(hdrs_dt, bins_dt, nomatch = 0)

# Count unique strains per bin
counts_PB <- overlaps[, .(n_strains = uniqueN(STRAIN)), by = .(CHROM, start, end)]

# If you want to merge with the full bin list (including 0s):
bins_wCounts <- merge(bins_dt, counts_PB, by = c("CHROM", "start", "end"), all.x = TRUE)
bins_wCounts[is.na(n_strains), n_strains := 0]

bins_wFreq <- as.data.frame(bins_wCounts) %>%
  dplyr::mutate(freq=n_strains/140) #change me to number of isotypes




# Merged SV plot
jasmine_plt <- ggplot(data = merged_SV %>% dplyr::filter(chrom != "MtDNA") %>% dplyr::select(chrom, pos, sv_type, number_svs_merged)) +
  geom_rect(aes(xmin = (pos / 1e6) - 0.001, xmax = (pos / 1e6) + 0.001, ymin = 0, ymax = number_svs_merged, fill = sv_type)) + 
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 14, color = 'black'),
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
    axis.title = element_text(size = 14, color = 'black'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0,0)) +
  labs(x = "Number of strains contributing to merged SV")
jasmine_plt_hist


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


snps_plt <- ggplot(snps) + 
  geom_point(aes(x = bin, y = variant_count), color = '#DB6333', alpha = 0.7) +
  facet_wrap( ~chrom, nrow = 1, scales = "free_x") + 
  # geom_smooth(aes(x = bin, y = variant_count), method = "loess", se = TRUE, color = "lightblue") +
  # ylab("Variants per kb") + 
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, color = 'black'),
    panel.grid = element_blank(),
    # axis.title.y = element_text(size = 14, color = 'black'),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size = 16, color = "black"))
snps_plt 


# # Funny supernova plot
# snps_circ <- snps %>%
#   arrange(chrom, bin) %>%
#   group_by(chrom) %>%
#   mutate(chr_len = max(bin)) %>%
#   ungroup() %>%
#   mutate(chr_offset = lag(cumsum(chr_len), default = 0)) %>%
#   mutate(x_circ = bin + chr_offset)
# 
# chr_labs <- snps_circ %>%
#   group_by(chrom) %>%
#   summarise(
#     center = mean(range(x_circ))
#   )
# 
# snps_plt <- ggplot(snps_circ) +
#   geom_point(
#     aes(x = x_circ, y = variant_count),
#     color = "#DB6333",
#     alpha = 0.7,
#     size = 0.6
#   ) +
#   coord_polar(theta = "x") +
#   scale_x_continuous(
#     breaks = chr_labs$center,
#     labels = chr_labs$chrom
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     axis.text.x = element_text(size = 14),
#     plot.margin = margin(10,10,10,10)
#   )
# snps_plt


