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
library(GenomicRanges)
library(IRanges)
library(ggrepel)

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
  dplyr::select(-IDY) %>% dplyr::filter(N2_chr == "IV") %>% dplyr::filter(strain != "ECA396") %>% 
  dplyr::filter(N2E >= 11000000) %>%
  dplyr::group_by(strain, contig) %>%
  dplyr::mutate(aln_span = sum(L2)) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(strain) %>%
  dplyr::filter(aln_span == max(aln_span)) %>%
  dplyr::ungroup() 

ecaAll <- ggplot(eca) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  geom_vline(xintercept = 13.46, color = 'black', linewidth = 2) +
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




# Looking at an INV that is found among 88 strains
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
merged_SV <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/pav/jasmine_SVmerging/elegans/output/clean_vcf.tsv")
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
geo_initial <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo.tsv")
hawaii_islands <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/elegans_isotypes_sampling_geo_hawaii_islands.tsv") %>% dplyr::select(isotype,collection_island_Hawaii)
WSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/140_correct.tsv", col_names = "strain") %>% dplyr::pull()
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
merged_all <- merged_SV %>% dplyr::filter(number_svs_merged == "141") %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::arrange(desc(sv_length)) %>%
  dplyr::ungroup() 
  # dplyr::slice_head(n = 3) %>%
  # dplyr::ungroup()
# 
# 
nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain != "ECA396")






# Look at SNVs inside and outside of INV for three Kaua'i strains!!
inv_rect <- filt_calls %>%
  dplyr::select(chrom,pos,sv_type, sv_length, strain) %>%
  dplyr::group_by(sv_type) %>%
  dplyr::arrange(desc(sv_length)) %>%
  dplyr::slice_head(n = 4) %>%
  dplyr::filter(sv_type == "INV") %>% dplyr::filter(strain != "JU3226")


inv <- nucmer %>% dplyr::filter(strain == "ECA1413" | strain == "ECA2968" | strain == "XZ1515") %>% dplyr::filter(N2_chr == "IV") #%>% dplyr::filter(contig == "ptg000006l")

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
  coord_cartesian(xlim = c(6.281940, 8)) +
  ggtitle("Largest inversion on CHROM IV") +
  labs(x = "N2 genome position (Mb)", y = "Wild strain contig position (Mb)")
inv1


genome_wide_haplotypeSharing <- nucmer %>% dplyr::filter(strain == "ECA1413" | strain == "ECA2968" | strain == "XZ1515")

ECA1413 <- ggplot(genome_wide_haplotypeSharing %>% dplyr::filter(strain == "ECA1413" & N2_chr != "MtDNA")) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
    facet_wrap(~N2_chr, nrow = 1, scales = "free") +
    theme(
      legend.position = 'none',
      axis.text = element_text(size = 14, color = 'black'),
      axis.ticks = element_blank(),
      axis.title.y = element_text(size = 16, color = 'black', face = 'bold'),
      axis.title.x = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA)) +
    labs(x = "N2 genome position (Mb)", y = "ECA1413")
ECA1413

ECA2968 <- ggplot(genome_wide_haplotypeSharing %>% dplyr::filter(strain == "ECA2968" & N2_chr != "MtDNA")) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  facet_wrap(~N2_chr, nrow = 1, scales = "free") +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)) +
  labs(x = "N2 genome position (Mb)", y = "ECA2968")
ECA2968

XZ1515 <- ggplot(genome_wide_haplotypeSharing %>% dplyr::filter(strain == "XZ1515" & N2_chr != "MtDNA")) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  facet_wrap(~N2_chr, nrow = 1, scales = "free") +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title.y = element_text(size = 16, color = 'black', face = 'bold'),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)) +
  labs(x = "N2 genome position (Mb)", y = "XZ1515")
XZ1515

kauai <- cowplot::plot_grid(
  ECA1413,XZ1515,ECA2968,
  nrow = 3,
  align = "v"
)
kauai


snvEC14 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/misc/PAV_INV_SNV_density/ECA1413_SNV_chromIV_INV.tsv", col_names = c("chrom",'pos','ref','alt')) %>%
  dplyr::mutate(strain = "ECA1413") %>%
  dplyr::filter(pos >= 6281940)
snvEC29 <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/misc/PAV_INV_SNV_density/ECA2968_SNV_chromIV_INV.tsv", col_names = c('chrom','pos','ref','alt')) %>%
  dplyr::mutate(strain = "ECA2968") %>%
  dplyr::filter(pos >= 6281940)
snvXZ <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/misc/PAV_INV_SNV_density/XZ1515_SNV_chromIV_INV.tsv", col_names = c("chrom",'pos','ref','alt')) %>%
  dplyr::mutate(strain = "XZ1515") %>%
  dplyr::filter(pos >= 6281940)
snvs_inINV <- snvEC14 %>% dplyr::bind_rows(snvEC29, snvXZ) %>% dplyr::mutate(pos = as.numeric(pos))

ggplot() + 
  # geom_point(data = snvs_inINV, aes(x = pos/1e6, y = 0.5), size = 2, color = 'black') +
  geom_density(data = snvs_inINV, aes(x = pos/1e6, alpha = 0.6), adjust = 0.5, fill = "firebrick") +
  geom_rect(data = inv_rect, aes(xmin = (pos + sv_length) / 1e6, xmax = pos / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
  facet_wrap(~strain, ncol = 1)+
  theme(
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(6.281940, 8)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.01)) +
  ggtitle("Largest inversion on CHROM IV") +
  labs(x = "N2 genome position (Mb)")



# SNVs inside and outside of massive INV in ECA3088
eca_all <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY) %>% dplyr::filter(strain == "ECA3088" & N2_chr == "IV" & contig == "ptg000002l")


ecaFour <- ggplot(eca_all) +
  geom_rect(aes(xmin = 13.46, xmax = 17492403 / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.01) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  facet_wrap(~N2_chr, scales = 'free') +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA)) +
  coord_cartesian(xlim = c(12, 18), ylim = c(8,22)) +
  labs(x = "N2 genome position (Mb)", y = "ECA3088 contig position (Mb)")
ecaFour

eca3088_snvs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/misc/PAV_INV_SNV_density/ECA3088_snvs_inINV.tsv",  col_names = c('chrom','pos','ref','alt', 'gt'))

ggplot() + 
  geom_rect(aes(xmin = 13.46, xmax = 17492403 / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.1) +
  # geom_point(data = eca3088_snvs, aes(x = pos/1e6, y = 0.2), size = 2, color = 'black') +
  geom_density(data = eca3088_snvs, aes(x = pos/1e6, alpha = 0.6), adjust = 0.5, fill = "firebrick") +
  geom_rect(aes(xmin = 13.5, xmax = 17492403 / 1e6, ymin = -Inf, ymax = Inf), fill = "gold", alpha = 0.5) +
  theme(
    axis.text = element_text(size = 14, color = 'black'),
    axis.ticks = element_blank(),
    legend.position = 'none',
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16, color = 'black', face = 'bold'),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 24, color = 'black', face = 'bold', hjust = 0.5)) +
  coord_cartesian(xlim = c(11, 18)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("ECA3088 inversion on CHROM IV") +
  labs(x = "N2 genome position (Mb)")










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


# Calculating HDR frequency
chr_order <- c("I","II","III","IV","V","X")
chrom_sizes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/N2.WS283.cleaned.fa.fai", col_names = c("chrom","start","end")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chr_order)) %>%
  dplyr::mutate(chrom = as.character(chrom))

bin_size <- 1000L

bins <- chrom_sizes %>%
  dplyr::group_by(chrom) %>%
  dplyr::do({
    chr_len <- .$end[1]
    starts <- seq(0, chr_len - 1, by = bin_size)
    tibble(
      chrom = .$chrom[1],
      start = starts,
      end   = pmin(starts + bin_size, chr_len)
    )
  }) %>%
  dplyr::ungroup()


# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
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
# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
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
# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
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
# bins <- snps %>% dplyr::select(chrom, bin)%>% dplyr::group_by(chrom) %>% dplyr::mutate(binEnd = lead(bin)) %>% dplyr::slice(-dplyr::n()) %>% dplyr::ungroup() %>% dplyr::rename(start = bin, end = binEnd)
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


# all_three <- cowplot::plot_grid(
#   del_freq, ins_freq, inv_freq,
#   nrow = 3,
#   align = "v",
#   rel_heights = c(1,1,1))
# all_three





# =======================================================#
# Making circlize plot #
# =======================================================#

# Outer ring of chromosomes (gene map of gene models represented with black rectangles)
## Chromosome IDs and sizes (start is always equal to zero)
chr_order <- c("I","II","III","IV","V","X")
chrom_sizes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/misc/N2.WS283.cleaned.fa.fai", col_names = c("chrom","start","end")) %>%
  dplyr::mutate(chrom = factor(chrom, levels = chr_order)) %>%
  dplyr::mutate(chrom = as.character(chrom))

## N2 gene modesl in BED format
n2_genes_bed <- n2_genes_plt %>% dplyr::filter(chrom != "MtDNA") %>%
  transmute(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end))

gene_bin_size <- 50000L

# one row per chr with lengths
chr_len_df <- chrom_sizes %>%
  dplyr::filter(chrom %in% chr_order) %>%
  dplyr::transmute(chrom = as.character(chrom), chr_len = as.numeric(end)) %>%
  dplyr::distinct()

# bins that cover the entire chromosome, including the last partial bin
gene_bins <- chr_len_df %>%
  dplyr::group_by(chrom) %>%
  dplyr::do({
    L <- .$chr_len[1]
    starts <- seq(0, L, by = gene_bin_size)  # include L so last bin is created
    tibble(
      chrom = .$chrom[1],
      start = starts[-length(starts)],
      end   = pmin(starts[-1], L)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(mid = (start + end)/2)

# overlap-count genes per bin
bins_dt  <- as.data.table(gene_bins)
genes_dt <- as.data.table(n2_genes_bed)  # chrom/start/end already numeric in your code

setkey(bins_dt,  chrom, start, end)
setkey(genes_dt, chrom, start, end)

ov <- foverlaps(genes_dt, bins_dt, nomatch = 0)

gene_counts <- ov[, .(gene_count = .N), by = .(chrom, start, end, mid)]

bins_gene <- merge(bins_dt, gene_counts, by = c("chrom","start","end","mid"), all.x = TRUE)
bins_gene[is.na(gene_count), gene_count := 0]

gene_bins_50kb <- as.data.frame(bins_gene) %>%
  dplyr::rename(pos = mid, value = gene_count)

gene_ylim <- c(0, max(gene_bins_50kb$value, na.rm = TRUE))








# Next ring in will have HDRs (represented with gray40 rectangles)
## HDRs in BED format
hdrs_collapsed_bed <- all_collapsed %>%
  transmute(chrom = as.character(chrom), start = as.numeric(start), end = as.numeric(end))

# Next ring in will have SNV count per kb (plot fitted LOESS line?) 
## SNV count in table: "chrom", "bin", "variant_count"
bin_size <- 1000L

bins <- chrom_sizes %>%
  dplyr::group_by(chrom) %>%
  dplyr::do({
    chr_len <- .$end[1]
    starts <- seq(0, chr_len - 1, by = bin_size)
    tibble(
      chrom = .$chrom[1],
      start = starts,
      end   = pmin(starts + bin_size, chr_len)
    )
  }) %>%
  dplyr::ungroup()

bins_dt <- as.data.table(bins)
setkey(bins_dt, chrom, start, end)
bins_dt[, id := .I]

snps_counts <- snps %>%
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(chrom %in% chr_order) %>%
  dplyr::mutate(
    chrom = as.character(chrom),
    start = (pos %/% bin_size) * bin_size,
    end   = start + bin_size) %>%
  dplyr::count(chrom, start, end, name = "variant_count")

snps_dt <- as.data.table(snps_counts)
setkey(snps_dt, chrom, start, end)

bins_snp <- merge(bins_dt, snps_dt, by = c("chrom","start","end"), all.x = TRUE)
bins_snp[is.na(variant_count), variant_count := 0]
snps_per_bin <- as.data.frame(bins_snp) %>%
  dplyr::mutate(pos = (start + end)/2) %>%
  dplyr::select(chrom,pos,variant_count) %>% dplyr::rename(value = variant_count)

# Next ring in will plot DEL frequncy (plotted as LOESS line) - make sure CGC1 is filtered out
## DEL calls with "chrom", "middle" (of the 1 kb bin), and "freq"
del_bin_freq <- del_bin_plt %>% dplyr::select(chrom, middle, freq)

# Next ring in will plot INS frequency (plotted as LOESS line) - make sure CGC1 is filtered out
## INS calls with "chrom", "middle" (of the 1 kb bin), and "freq"
ins_bin_freq <- ins_bin_plt %>% dplyr::select(chrom, middle, freq)

# Then the final, more inner ring will have INV frequency (plotted as LOESS line) - make sure CGC1 is filtered out
## INV calls with "chrom", "middle" (of the 1 kb bin), and "freq"
inv_bin_freq <- inv_bin_plt %>% dplyr::select(chrom, middle, freq)


## =========================
##  Plotting!
## ========================= 
# SV frequency tracks must have: chrom, pos, value
del_bin_freq <- del_bin_freq %>%
  transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))

ins_bin_freq <- ins_bin_freq %>%
  transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))

inv_bin_freq <- inv_bin_freq %>%
  transmute(chrom = as.character(chrom), pos = as.numeric(middle), value = as.numeric(freq))

add_rect_track <- function(bed, col, track_height = 0.08, label = NULL) {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- bed[bed$chrom == chr, , drop = FALSE]
      if (nrow(d) == 0) return()
      circos.rect(d$start, 0, d$end, 1, col = col, border = NA)
    }
  )
  
  if (!is.null(label)) {
    circos.text(
      x = 0, y = 0.5, labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}

add_value_track <- function(df, col, ylim, track_height = 0.10, label = NULL,
                            type = c("line","points"), add_loess = FALSE, span = 0.2) {
  type <- match.arg(type)
  
  circos.trackPlotRegion(
    ylim = ylim,
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- df[df$chrom == chr, , drop = FALSE]
      if (nrow(d) == 0) return()
      
      d <- d[order(d$pos), , drop = FALSE]
      
      if (type == "points") {
        circos.points(d$pos, d$value, pch = 16, cex = 0.25, col = col)
      } else {
        circos.lines(d$pos, d$value, col = col, lwd = 1)
      }
      
      if (add_loess && nrow(d) >= 50) {
        fit <- stats::loess(value ~ pos, data = d, span = span)
        xs  <- d$pos
        ys  <- stats::predict(fit, newdata = data.frame(pos = xs))
        ok  <- is.finite(ys)
        if (any(ok)) circos.lines(xs[ok], ys[ok], col = col, lwd = 2)
      }
    }
  )
  
  if (!is.null(label)) {
    circos.text(
      x = 0, y = mean(ylim), labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}

add_filled_area_track <- function(df, col_fill = "#00BFC480", col_line = "#00BFC4",
                                  ylim, track_height = 0.10, label = NULL) {
  circos.trackPlotRegion(
    ylim = ylim,
    track.height = track_height,
    bg.border = NA,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      d <- df[df$chrom == chr, , drop = FALSE]
      if (nrow(d) < 2) return()
      d <- d[order(d$pos), , drop = FALSE]
      
      # build polygon (baseline at ylim[1])
      x_poly <- c(d$pos, rev(d$pos))
      y_poly <- c(d$value, rep(ylim[1], nrow(d)))
      
      circos.polygon(x_poly, y_poly, col = col_fill, border = NA)
      circos.lines(d$pos, d$value, col = col_line, lwd = 1.2)
    }
  )
  
  if (!is.null(label)) {
    circos.text(
      x = 0, y = mean(ylim), labels = label,
      sector.index = chr_order[1],
      track.index  = get.current.track.index(),
      facing = "inside", adj = c(1, 0.5), cex = 0.7
    )
  }
}

circos.clear()
circos.par(track.margin = c(0.002, 0.002), cell.padding = c(0, 0, 0, 0)) # start.degree = 86, gap.after = c(rep(2, length(chr_order)-1), 8)
circos.initialize(factors = as.character(chrom_sizes$chrom), xlim = cbind(rep(0, nrow(chrom_sizes)), chrom_sizes$end))

# Outer ideogram-like ring + gene models inside it
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.02,
  bg.border = NA,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    
    # background chromosome band
    circos.rect(xlim[1], 0, xlim[2], 1, col = "grey90", border = "black")
    # gene models as black rectangles (thin band)
    # g <- n2_genes_bed[n2_genes_bed$chrom == chr, , drop = FALSE]
    # if (nrow(g) > 0) {
    #   circos.rect(g$start, 0, g$end, 1, col = "black", border = NA)}
    
    # gene density (50 kb bins) as cyan-blue line
    # gd <- gene_density_50kb[gene_density_50kb$chrom == chr, , drop = FALSE]
    # if (nrow(gd) > 1) {
    #   gd <- gd[order(gd$pos), , drop = FALSE]
    #   
    #   # rescale density to fit nicely within the outer track band
    #   y <- (gd$value - gene_ylim[1]) / (gene_ylim[2] - gene_ylim[1] + 1e-9)
    #   y <- 0.10 + y * 0.80   # occupy 10%..90% of the track height
    #   
    #   circos.lines(gd$pos, y, col = "#00BFC4", lwd = 2)
    # }
    
    
    
    # chromosome label
    circos.text(CELL_META$xcenter, 3.3, chr, facing = "bending.outside", niceFacing = T, cex = 1.5, font = 2)
    # axis ticks
    circos.axis(h = "top", major.at = seq(0, xlim[2], by = 5e6),labels = seq(0, xlim[2], by = 5e6) / 1e6, labels.cex = 1, major.tick.length = 0.001)
  }
)


add_filled_area_track(gene_bins_50kb,col_fill = "#00BFC480",col_line = "#00BFC4", ylim = gene_ylim,track_height = 0.10)
add_rect_track(hdrs_collapsed_bed, col = "grey40", track_height = 0.08)
snp_ylim <- c(0, max(snps_per_bin$value, na.rm = TRUE))
add_value_track(snps_per_bin, col = "#DB6333", ylim = snp_ylim, track_height = 0.12, type = "line", add_loess = FALSE, span = 0.15)
add_value_track(del_bin_freq, col = "red",  ylim = c(0, 1), track_height = 0.08, type = "line", add_loess = FALSE, span = 0.2)
add_value_track(ins_bin_freq, col = "blue", ylim = c(0, 1), track_height = 0.08, type = "line", add_loess = FALSE, span = 0.2)
add_value_track(inv_bin_freq, col = "gold3", ylim = c(0, 1), track_height = 0.08, type = "line", add_loess = FALSE, span = 0.2)

lgd <- Legend(title = "Tracks",labels = c("Genes per 50 kb","HDRs","SNPs per kb","DEL frequency","INS frequency","INV frequency"),
              type = "lines",legend_gp = gpar(col = c("#00BFC4","grey40","#DB6333","red","blue","gold3"),
              lwd = c(3, 3, 3, 3, 3, 3)),
              labels_gp = gpar(fontsize = 14),
              title_gp  = gpar(fontsize = 16, fontface = "bold"))

draw(lgd, x = unit(0.92, "npc"), y = unit(0.92, "npc"), just = c("right", "top"))

circos.clear()




























# =============================================================== #
# Enrichment of SVs in HDRs #
# =============================================================== #
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

# all_SVs <- merged_SV %>% 
#   dplyr::select(chrom,pos,sv_length,sv_type,number_svs_merged) %>%
#   dplyr::mutate(sv_length = abs(sv_length)) %>%
#   dplyr::mutate(end = pos + sv_length) %>% 
#   dplyr::rename(start = pos) %>%
#   dplyr::select(chrom, start, end, sv_type) 

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



################## GPCRS! #####################
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





################## Histones.... #####################
ipr_test <- ipr %>% dplyr::select(IPR_description) %>% dplyr::distinct()
# ipr_histone <- ipr %>% dplyr::filter(grepl("Histone", IPR_description) & !grepl('ase', IPR_description) & !grepl('link', IPR_description) & !grepl('fold', IPR_description)) %>% dplyr::select(N2, IPR_description) %>% dplyr::distinct(N2, .keep_all = T)
ipr_histone <- ipr %>% dplyr::filter(grepl("Histone H", IPR_description) & !grepl('ase', IPR_description) & !grepl("CENP", IPR_description)) %>% dplyr::select(N2, IPR_description) %>% dplyr::distinct(N2, .keep_all = T)
# ipr_hist_terms <- ipr_histone %>% dplyr::distinct(IPR_description)

histones <- cleaned_N2_gff %>%
  dplyr::left_join(ipr_histone, by = "N2") %>% 
  dplyr::select(chrom, start, end, N2, IPR_description) %>% 
  dplyr::filter(!is.na(IPR_description))

histone_list <- histones %>% dplyr::pull(N2)
total_histone <- length(histone_list)

# gpc_test <- ipr %>% dplyr::filter(grepl("7TM", IPR_description)) %>%
#   dplyr::group_by(N2) %>%
#   dplyr::mutate(annotation_count = n()) %>%
#   dplyr::ungroup()
#  
# multiple <- gpc_test %>% dplyr::filter(annotation_count > 1) %>% dplyr::distinct(N2)
#  
# at_least_one_anno <- gpc_test %>% dplyr::distinct(N2)

svs_dt <- as.data.table(common_SVs)
histone_dt <- as.data.table(histones)

setkey(svs_dt, chrom, start, end)
setkey(histone_dt, chrom, start, end)

svs_inCodingRegions <- data.table::foverlaps(x = svs_dt, y = histone_dt, type = "any") %>% dplyr::filter(!is.na(start)) %>% dplyr::mutate(overlap = T)

overlap <- svs_inCodingRegions %>% dplyr::select(chrom, start, end, N2, sv_type, overlap) %>% dplyr::distinct(chrom,start,end,N2,sv_type, .keep_all = T)

plt_stats_juvenile <- overlap %>%
  dplyr::filter(N2 %in% histone_list) %>%
  dplyr::distinct(N2, sv_type) %>%          
  dplyr::count(sv_type, name = "n_genes") %>%
  dplyr::mutate(total_histone = total_histone, proportion = 100 * n_genes / total_histone) %>%
  dplyr::mutate(region = T)

plt_stats <- plt_stats_juvenile %>%
  tidyr::complete(sv_type, fill = list(n_genes = 0, region = "TRUE", total_histone = total_histone)) %>%
  dplyr::mutate(n_false = total_histone - n_genes) %>%
  select(sv_type, n_true = n_genes, n_false, total_histone) %>%
  pivot_longer(c(n_true, n_false), names_to = "region", values_to = "n_genes") %>%
  dplyr::mutate(region = recode(region, n_true = "TRUE", n_false = "FALSE"),proportion = 100 * n_genes / total_histone)

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



















# ============================================== # 
# SV's correlation with geography
# ============================================== # 
MAF_thresh <- round(0.05 * 141)

SVs_geo <- merged_SV %>% 
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(number_svs_merged > MAF_thresh) %>%
  dplyr::mutate(sv_length = abs(sv_length)) %>%
  dplyr::mutate(end = pos + sv_length) %>% 
  dplyr::rename(start = pos) %>%
  tidyr::pivot_longer(
       cols = -c(chrom, start, end, sv_type, sv_length, number_svs_merged),  # Keep chrom, pos, and info as identifiers
       names_to = "strain",          # New column for sample names
       values_to = "genotype"           # New column for sample values
       ) %>%
     dplyr::mutate(genotype=ifelse(genotype=="./.",0,1)) %>%
  dplyr::left_join(geo, by = c("strain" = "isotype")) %>%
  dplyr::select(chrom, start, end, sv_type, sv_length, number_svs_merged, strain, genotype, geo) %>%
  dplyr::filter(genotype != "0") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "England", geo)) %>%
  dplyr::mutate(geo2 = ifelse(geo == "Oahu" | geo == "Big Island" | geo == "Maui" | geo == "Kauai", "Hawaii", geo)) 

############### NO ISLAND RESOLUTION ON HAWAII ######################################
clustering <- SVs_geo %>%
  dplyr::group_by(chrom,start,end,geo2) %>%
  dplyr::mutate(geo_max_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom,start,end) %>%
  dplyr::mutate(prop = max(geo_max_count) / dplyr::n()) %>%
  dplyr::ungroup()

cleaned <- clustering %>% 
  dplyr::select(chrom, start, sv_type, sv_length, number_svs_merged, geo, geo2, geo_max_count, prop) %>%
  dplyr::group_by(chrom, start, sv_length) %>%
  dplyr::filter(geo_max_count == max(geo_max_count)) %>%
  dplyr::distinct(chrom, start, sv_length,geo2, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom,start,sv_length) %>%
  dplyr::mutate(number_in_group = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(final_geo = ifelse(number_in_group > 1, "max_in_multiple_locations",geo2))

final <- cleaned %>% 
  dplyr::select(chrom, start, sv_type, sv_length, number_svs_merged, geo_max_count, prop, final_geo) %>%
  dplyr::distinct() %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS", "DEL", "INV"))) %>%
  dplyr::mutate(final_geo = factor(final_geo, levels = c("Africa","Atlantic","Europe","Hawaii","North America","Oceania","max_in_multiple_locations")))


geo.colors <- c("Hawaii"="#66C2A5", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="orange", "max_in_multiple_locations" = 'grey')

sv_corr_geo <- ggplot(data = final) +
  geom_point(aes(x = start / 1e6, y = prop, fill = final_geo, shape = sv_type), size = 2) + 
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
    ) +
  labs(x = "N2 genome coordinates (Mb)", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))),fill = guide_legend(override.aes = list(shape = 21, color = 'black', size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01))
sv_corr_geo

sv_corr_geo_noH <- ggplot(data = final %>% dplyr::filter(final_geo != "Hawaii")) +
  geom_point(aes(x = start / 1e6, y = prop, fill = final_geo, shape = sv_type), size = 2) + 
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "N2 genome coordinates (Mb)", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))),fill = guide_legend(override.aes = list(shape = 21, color = 'black',size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01))
sv_corr_geo_noH


final_freq <- final %>% dplyr::mutate(MAF = number_svs_merged / 141)
  
sv_corr_freq <- ggplot(data = final_freq) +
  geom_point(aes(x = MAF, y = prop, fill = final_geo, shape = sv_type), size = 2) +
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "Allele frequency", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))), fill = guide_legend(override.aes = list(shape = 21, color = 'black', size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01))
sv_corr_freq


norm_freq <- SVs_geo %>% dplyr::select(strain,geo2) %>% dplyr::distinct() %>% dplyr::group_by(geo2) %>% dplyr::mutate(geo_count = n()) %>% dplyr::ungroup()
final_norm_freq <- final_freq %>% dplyr::left_join(norm_freq, by = c("final_geo"= "geo2")) %>% dplyr::mutate(norm_af = geo_max_count / geo_count) %>%
  dplyr::filter(final_geo != "max_in_multiple_locations")

sv_corr_freq_norm <- ggplot(data = final_norm_freq) +
  geom_point(aes(x = norm_af, y = prop, fill = final_geo, shape = sv_type), size = 2) +
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +s
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "Normalized allele frequency", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))), fill = guide_legend(override.aes = list(shape = 21, color = 'black', size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01))
sv_corr_freq_norm  


############################ ISLAND RESOLUTION ON HAWAII ##############################################
clustering <- SVs_geo %>%
  dplyr::select(-geo2) %>%
  dplyr::group_by(chrom,start,end,geo) %>%
  dplyr::mutate(geo_max_count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom,start,end) %>%
  dplyr::mutate(prop = max(geo_max_count) / dplyr::n()) %>%
  dplyr::ungroup()

cleaned <- clustering %>% 
  dplyr::select(chrom, start, sv_type, sv_length, number_svs_merged, geo, geo_max_count, prop) %>%
  dplyr::group_by(chrom, start, sv_length) %>%
  dplyr::filter(geo_max_count == max(geo_max_count)) %>%
  dplyr::distinct(chrom, start, sv_length,geo, .keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(chrom,start,sv_length) %>%
  dplyr::mutate(number_in_group = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(final_geo = ifelse(number_in_group > 1, "max_in_multiple_locations",geo))

final <- cleaned %>% 
  dplyr::select(chrom, start, sv_type, sv_length, number_svs_merged, geo_max_count, prop, final_geo) %>%
  dplyr::distinct() %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(sv_type = factor(sv_type, levels = c("INS", "DEL", "INV"))) %>%
  dplyr::mutate(final_geo = factor(final_geo, levels = c("Africa","Atlantic","Big Island","Europe","Kauai","Maui", "Molokai", "North America","Oahu","Oceania","max_in_multiple_locations")))


geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="orange", "max_in_multiple_locations" = 'grey')

sv_corr_geo <- ggplot(data = final) +
  geom_point(aes(x = start / 1e6, y = prop, fill = final_geo, shape = sv_type), size = 2) + 
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "N2 genome coordinates (Mb)", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))),fill = guide_legend(override.aes = list(shape = 21, color = 'black', size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01))
sv_corr_geo

sv_corr_geo_noH <- ggplot(data = final %>% dplyr::filter(final_geo != "Oahu" & final_geo != "Maui" & final_geo != "Big Island" & final_geo != "Kauai")) +
  geom_point(aes(x = start / 1e6, y = prop, fill = final_geo, shape = sv_type), size = 2) + 
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "N2 genome coordinates (Mb)", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))),fill = guide_legend(override.aes = list(shape = 21, color = 'black',size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01))
sv_corr_geo_noH


final_freq <- final %>% dplyr::mutate(MAF = number_svs_merged / 141)

sv_corr_freq <- ggplot(data = final_freq) +
  geom_point(aes(x = MAF, y = prop, fill = final_geo, shape = sv_type), size = 2) +
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "Allele frequency", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))), fill = guide_legend(override.aes = list(shape = 21, color = 'black', size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01))
sv_corr_freq



norm_freq <- SVs_geo %>% dplyr::select(strain,geo) %>% dplyr::distinct() %>% dplyr::group_by(geo) %>% dplyr::mutate(geo_count = n()) %>% dplyr::ungroup()
final_norm_freq <- final_freq %>% dplyr::left_join(norm_freq, by = c("final_geo"= "geo")) %>% dplyr::mutate(norm_af = geo_max_count / geo_count) %>%
  dplyr::filter(final_geo != "max_in_multiple_locations")

sv_corr_freq_norm <- ggplot(data = final_norm_freq) +
  geom_point(aes(x = norm_af, y = prop, fill = final_geo, shape = sv_type), size = 2) +
  scale_fill_manual(values = geo.colors) +
  scale_shape_manual(values = c("INS" = 22, "DEL" = 23, "INV" = 24)) +
  facet_wrap(~chrom, scales = "free_x") +
  theme(
    # panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 13, color = 'black'),
  ) +
  labs(x = "Normalized allele frequency", y = "Proportion of strain geographic isolation with ALT allele / all strains with ALT allele", fill = "Isolation location", size = "SV type") +
  guides(size = guide_legend(override.aes = list(size = c(2, 4, 6))), fill = guide_legend(override.aes = list(shape = 21, color = 'black', size = 7))) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01))
sv_corr_freq_norm 









# ========================================================================================= #
# PCA on SVs
# ========================================================================================= #
# Filter for MAF > 0.05
n <- 141
maf <- 0.05
least <- ceiling(maf * n)      #8
most  <- floor((1 - maf) * n)  #133

common_vcf <- merged_SV %>% 
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most) %>%
  dplyr::select(-chrom, -pos, -ref, -alt, -sv_type, -sv_length, -number_svs_merged) %>%
  as.matrix()

common_vcf[common_vcf == "./."] <- 0
sv_mat <- apply(common_vcf, 2, as.numeric)

sv_mat_t <- t(sv_mat)

sv_mat_t[is.na(sv_mat_t)] <- colMeans(sv_mat_t, na.rm = TRUE)

sv_pca <- prcomp(sv_mat_t, center = TRUE, scale. = TRUE) 


pca_df <- as.data.frame(sv_pca$x)
pca_df$strain <- rownames(pca_df)

strain_geo <- geo %>% dplyr::rename(strain = isotype) %>% dplyr::select(strain,geo)

pca_df <- pca_df %>%
  dplyr::left_join(strain_geo, by = "strain") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "CGC1",geo)) 

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "All SVs", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 


dim(sv_mat_t)
summary(rowSums(sv_mat_t))
cor(rowSums(sv_mat_t), sv_pca$x[,1])  # PC1 is 93% correlated with number of SVs per strain


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "All SVs", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "black") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")


### DELETIONS ###
common_vcf <- merged_SV %>% 
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most & sv_type == "DEL") %>%
  dplyr::select(-chrom, -pos, -ref, -alt, -sv_type, -sv_length, -number_svs_merged) %>%
  as.matrix()

common_vcf[common_vcf == "./."] <- 0
sv_mat <- apply(common_vcf, 2, as.numeric)

sv_mat_t <- t(sv_mat)

sv_mat_t[is.na(sv_mat_t)] <- colMeans(sv_mat_t, na.rm = TRUE)

sv_pca <- prcomp(sv_mat_t, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(sv_pca$x)
pca_df$strain <- rownames(pca_df)

strain_geo <- geo %>% dplyr::rename(strain = isotype) %>% dplyr::select(strain,geo)

pca_df <- pca_df %>%
  dplyr::left_join(strain_geo, by = "strain") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "CGC1",geo)) 

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Deletions", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Deletions", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "red") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")





### INSERTIONS ###
common_vcf <- merged_SV %>% 
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most & sv_type == "INS") %>%
  dplyr::select(-chrom, -pos, -sv_type, -sv_length, -number_svs_merged) %>%
  as.matrix()

common_vcf[common_vcf == "./."] <- 0
sv_mat <- apply(common_vcf, 2, as.numeric)

sv_mat_t <- t(sv_mat)

sv_mat_t[is.na(sv_mat_t)] <- colMeans(sv_mat_t, na.rm = TRUE)

sv_pca <- prcomp(sv_mat_t, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(sv_pca$x)
pca_df$strain <- rownames(pca_df)

strain_geo <- geo %>% dplyr::rename(strain = isotype) %>% dplyr::select(strain,geo)

pca_df <- pca_df %>%
  dplyr::left_join(strain_geo, by = "strain") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "CGC1",geo)) 

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 50, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location", title = "Insertions", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 


# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Insertions", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "blue") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,17)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")



### INVERSIONS ###
common_vcf <- merged_SV %>% 
  dplyr::select(-ref,-alt) %>%
  dplyr::filter(number_svs_merged >= least & number_svs_merged <= most & sv_type == "INV") %>%
  dplyr::select(-chrom, -pos, -sv_type, -sv_length, -number_svs_merged) %>%
  as.matrix()

common_vcf[common_vcf == "./."] <- 0
sv_mat <- apply(common_vcf, 2, as.numeric)

sv_mat_t <- t(sv_mat)

sv_mat_t[is.na(sv_mat_t)] <- colMeans(sv_mat_t, na.rm = TRUE)

sv_pca <- prcomp(sv_mat_t, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(sv_pca$x)
pca_df$strain <- rownames(pca_df)

strain_geo <- geo %>% dplyr::rename(strain = isotype) %>% dplyr::select(strain,geo)

pca_df <- pca_df %>%
  dplyr::left_join(strain_geo, by = "strain") %>%
  dplyr::mutate(geo = ifelse(strain == "CGC1", "CGC1",geo)) 

geo.colors <- c("Big Island"="#66C2A5", "Molokai" = "black", "Maui" = "yellow", "Oahu" = "brown", "Kauai" = "pink", "Africa"="green", "North America" = "purple", "Europe" = "#E41A1C", "Atlantic" = "blue", 
                "Oceania" ="cyan", "unknown" = 'gray', "CGC1" = "#DB6333")

pca_df <- pca_df %>%
  mutate(label = ifelse(PC2 > 5, strain, NA))

ggplot(pca_df, aes(PC1, PC2, color = geo)) +
  geom_text_repel(aes(label = label), size = 4, max.overlaps = Inf, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Inversions", x = paste0("PC1 (", round(100 * summary(sv_pca)$importance[2,1], 1), "%)"),y = paste0("PC2 (", round(100 * summary(sv_pca)$importance[2,2], 1), "%)")) +
  guides(color = guide_legend(override.aes = list(size = 7))) 
  

# PC3/4
ggplot(pca_df, aes(PC3, PC4, color = geo)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = geo.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 12, color = 'black', face = 'bold'),
    legend.text = element_text(size = 14, color = 'black'),
    legend.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  labs(color = "Collection location",title = "Inversions", x = paste0("PC3 (", round(100 * summary(sv_pca)$importance[2,3], 1), "%)"),y = paste0("PC4 (", round(100 * summary(sv_pca)$importance[2,4], 1), "%)"))+
  guides(color = guide_legend(override.aes = list(size = 7))) 

# Scree plot:
scree <- data.frame(PC = seq_along(sv_pca$sdev),variance = sv_pca$sdev^2 / sum(sv_pca$sdev^2))

ggplot(scree[1:10,], aes(PC, variance * 100)) +
  geom_col(fill = "gold") +
  theme(
    axis.text = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    plot.margin = margin(10,10,10,10),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.title = element_text(size = 16, color = 'black', face = 'bold')) +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(breaks = 1:10) +
  coord_cartesian(ylim = c(0,16)) +
  labs(y = "Proportion of variance explained (%)", x = "Principal component")



