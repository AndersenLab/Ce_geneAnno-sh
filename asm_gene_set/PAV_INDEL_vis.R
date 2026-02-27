library(plyr) # ALWAYS LOAD BEFORE DPLYR
library(dplyr)
library(ggplot2)
library(readr)
library(IRanges)
library(tidyr)
library(valr)
library(data.table)


vcf <- readr::read_tsv("/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/pav_vcfs/pav_indels_invs_merged.tsv") %>%
  dplyr::rename(CHROM=`#CHROM`) %>%
  dplyr::filter(CHROM != "MtDNA")


ggplot() + 
  # geom_rect(data = all_collapsed %>% dplyr::rename(CHROM = chrom), aes(xmin = start/1e6, xmax = end/1e6, ymin = 0, ymax = 1.75, fill = "hdr")) + 
  geom_rect(data= vcf %>% dplyr::filter(SV == "DEL"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 1, ymax = 1.5, fill = "DEL")) +
  geom_rect(data = vcf %>% dplyr::filter(SV == "INS"), aes(xmin = POS/1e6, xmax = (POS + SV_length) / 1e6, ymin = 0.25, ymax = 0.75, fill = "INS")) +
  geom_rect(data = vcf %>% dplyr::filter(SV == "INV"), aes(xmin = POS/1e6, xmax = (POS + SV_length) / 1e6, ymin = 1.75, ymax = 2.25, fill = "INV")) +
  # geom_rect(data = )
  facet_wrap(~CHROM, nrow=6, scales = "free_x") + 
  scale_fill_manual(
    name = "",
    values = c("DEL" = "red", 'INS' = 'blue', 'INV' = 'gold'),
    breaks = c("Deletions", "Insertion")
  ) + 
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text("Genome position (Mb)", size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line.x = element_line(),
    legend.position = "none",
    strip.text = element_text(size = 14, color = 'black'),
    panel.background = element_blank(),
    # plot.margin = margin(0, 0, 0, 0)
  )

vcf_longer <- vcf %>%
  dplyr::select(-REF,-ALT) %>%
  tidyr::pivot_longer(
    cols = -c(CHROM, POS, SV, SV_length),  
    names_to = "sample",          
    values_to = "genotype"        
  ) %>%
  dplyr::mutate(genotype=ifelse(genotype=="./.",0,1)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n_alt=sum(genotype)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_alt)) %>%
  dplyr::mutate(rid=rleid(n_alt)) %>%
  dplyr::filter(genotype > 0) 

### Then you will need to intersect with HDRs for a particular strain and add a column to indicate if that variant is found in a HDR for that particular strain
SV_plt <- ggplot(vcf_longer) +
  geom_rect(data = vcf_longer %>% dplyr::filter(SV == "DEL"), aes(xmin=(POS-500)/1e6, xmax=(POS + 500)/1e6, ymin=rid+0.5, ymax=rid-0.5, fill=SV)) +
  geom_rect(data = vcf_longer %>% dplyr::filter(SV != "DEL"), aes(xmin=(POS-500)/1e6, xmax=(POS + SV_length + 500)/1e6, ymin=rid+0.5, ymax=rid-0.5, fill=SV)) + 
  facet_wrap(~CHROM,nrow=1,scales = 'free_x') +
  scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 13, color = 'black', face = 'bold'),
    legend.position = "none") +
  xlab("N2 Genome position (Mb)") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  ylab("115 wild strains")
SV_plt

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/PAV_allStrains.png", SV_plt, width = 15, height = 12, dpi = 600)



strain <- vcf_longer %>% dplyr::select(sample) %>% dplyr::distinct(sample) %>% dplyr::pull()
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom","start","end","sample")) %>%
  dplyr::arrange(sample) %>% 
  dplyr::filter(sample %in% strain)


# SV_plt_noHDR <- ggplot(vcf_longer) + 
#   geom_rect(data = hdrs %>% dplyr::filter(chrom == "V"), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.2, ymax = 0.8), fill = 'gray70') +
#   geom_rect(data = vcf_longer %>% dplyr::filter(CHROM == "V" & SV != "DEL"), aes(xmin = (POS - 500)/1e6, xmax = (POS + 500)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = vcf_longer %>% dplyr::filter(CHROM == "V" & SV == "DEL"), aes(xmin = (POS - 500)/1e6, xmax = (POS + 500)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# SV_plt_noHDR


ROI <- vcf_longer %>% dplyr::filter(CHROM == "V" & POS >= 19200000 & POS <= 19300000) %>% dplyr::rename(chrom = CHROM)
ROI_hdrs <- hdrs %>% 
  dplyr::filter(chrom == "V") 

# Convert to data.tables
ROI_dt <- as.data.table(ROI)
ROI_hdrs_dt <- as.data.table(ROI_hdrs)

# Create start and end columns for ROI
ROI_dt[, `:=`(start = POS, end = POS)]

# We want to match by sample AND overlap on positions
# So split by sample first (or use a join on sample)
merged <- ROI_hdrs_dt[ROI_dt, on = .(sample), allow.cartesian = TRUE]

# Now perform overlap check
plot_data <- merged[is.na(start) | (start <= POS & end >= POS)]

plt_data <- plot_data %>%
  dplyr::select(-i.chrom, -i.start, -i.end) %>%
  dplyr::rename(hdr_start = start, hdr_end = end) %>%
  dplyr::mutate(
    hdr_start_updated = ifelse(.data$hdr_start < 19200000, 19200000, .data$hdr_start),
    hdr_end_updated   = ifelse(.data$hdr_end > 19300000, 19300000, .data$hdr_end)) %>%
  dplyr::select(chrom,POS,SV,SV_length,sample,genotype,n_alt,rid,hdr_start_updated,hdr_end_updated)

plot_data_final <- ROI %>%
  dplyr::left_join(plt_data, by = c("POS","sample")) %>%
  dplyr::rename(chrom = chrom.x, SV = SV.x, n_alt = n_alt.x, SV_length = SV_length.x) %>%
  dplyr::select(chrom,POS,SV,SV_length,sample,hdr_start_updated,hdr_end_updated)


# Adding in rows for strains that have HDRs in the ROI, but so SVs in the ROI so the HDR calls are not being pulled
SV_samples <- unique(c(plot_data_final$sample))
missing_samples <- setdiff(strain, SV_samples)

ROIs <- 19200000
ROIe <- 19300000
ROI_hdrs_missing <- ROI_hdrs %>%
  dplyr::filter(sample %in% missing_samples) %>%
  dplyr::filter(
    start < ROIs & end > ROIe | 
      start < ROIs & end > ROIs |
      start < ROIe & end > ROIe) %>%
  dplyr::mutate(chrom ="V", POS = NA, SV = NA, SV_length = NA) %>%
  dplyr::rename(hdr_start = start, hdr_end = end) %>%
  dplyr::select(chrom, POS, SV, SV_length, sample, hdr_start, hdr_end) %>%
  dplyr::mutate(
    hdr_start_updated = ifelse(hdr_start < 19200000, 19200000, hdr_start),
    hdr_end_updated   = ifelse(hdr_end > 19300000, 19300000, hdr_end)) %>%
  dplyr::select(-hdr_start, -hdr_end)

# plt_data_full <- bind_rows(plot_data_final, ROI_hdrs_missing)




### Addingh in HDRs that are in ROI but not directly overlappign with a SV so they are not pulled in the initial merge step
hdr_dt <- as.data.table(ROI_hdrs)
var_dt <- as.data.table(plot_data_final %>% filter(!is.na(POS)))

# Add start and end columns for variants (positions as single-point ranges)
var_dt[, `:=`(var_start = POS, var_end = POS)]

# Key by sample and ranges
setkey(hdr_dt, sample, start, end)
setkey(var_dt, sample, var_start, var_end)

# Find overlaps
overlaps <- foverlaps(var_dt, hdr_dt, 
                      by.x = c("sample","var_start","var_end"), 
                      by.y = c("sample","start","end"),
                      type = "any", nomatch = 0L)

# These HDRs have at least one overlapping variant
hdr_with_variants[, overlap := TRUE]

# Merge to flag which HDRs have overlaps
HDR_no_overlap <- merge(
  hdr_dt, hdr_with_variants,
  by = c("sample","start","end"),
  all.x = TRUE
)[is.na(overlap)]  # keep only those with no overlap

# Filter to ROI and reformat for plotting
HDR_no_overlap <- HDR_no_overlap %>%
  dplyr::filter(chrom == "V", start < ROIe & end > ROIs) %>%
  dplyr::mutate(
    chrom = "V",
    POS = as.numeric(NA),
    SV = as.character(NA),
    SV_length = as.numeric(NA),
    hdr_start_updated = ifelse(start < ROIs, ROIs, start),
    hdr_end_updated = ifelse(end > ROIe, ROIe, end)
  ) %>%
  dplyr::select(chrom, POS, SV, SV_length, sample, hdr_start_updated, hdr_end_updated)

# Combine with the previous set
plt_data_full <- bind_rows(plot_data_final, ROI_hdrs_missing, HDR_no_overlap)



SV_plot <- ggplot(plt_data_full) + 
  geom_rect(data = plt_data_full %>% dplyr::filter(!is.na(hdr_start_updated)), aes(xmin = hdr_start_updated / 1e6, xmax = hdr_end_updated / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
  geom_rect(data = plt_data_full %>% dplyr::filter(!is.na(SV)) %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
  geom_rect(data = plt_data_full %>% dplyr::filter(!is.na(SV)) %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
  facet_wrap(~sample) +
  scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  scale_x_continuous(expand = c(0,0))
SV_plot

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/SV_ROI_stefan.png", SV_plot, width = 25, height = 15, dpi = 600)


# 
# 
# # Looking at specific ROI #
# # kola genes for robert: 
# genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
# all_genes_strain <- genes_strain %>%
#   dplyr::filter(strain != "N2" | grepl("protein_coding", attributes)) %>% 
#   dplyr::mutate(attributes = gsub("ID=gene:","", attributes)) %>%
#   dplyr::mutate(attributes = gsub("ID=","", attributes)) %>%
#   dplyr::mutate(attributes = sub(";.*", "", attributes)) 
# 
# n2_genes <- all_genes_strain %>% dplyr::filter(strain == "N2") %>% 
#   dplyr::filter(attributes %in% c("WBGene00001626", "WBGene00019435", "WBGene00019957", "WBGene00015280", "WBGene00008194", "WBGene00010605")) %>% 
#   dplyr::rename(chrom = contig)
# 
# 
# # K06A9.1 - WBGene00019435
# # X 1547505 - 1560145
# lowerk = 1525000
# upperk = 1580000
# K1 <- vcf_longer %>% dplyr::filter(CHROM == "X" & POS > lowerk & POS < upperk)
# K1_hdr <- hdrs %>% 
#   dplyr::filter(chrom == "X") %>%
#   dplyr::filter(start > lowerk & start < upperk & end > lowerk & end < upperk)
# 
# K1_SV_plt <- ggplot(K1) + 
#   geom_rect(data = K1_hdr, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
#   annotate("rect", xmin = 1547505 / 1e6, xmax = 1560145 / 1e6, ymin = -Inf, ymax = Inf, fill = 'seagreen', alpha = 0.3) +
#   geom_rect(data = K1 %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = K1 %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# K1_SV_plt
# 
# 
# # R08E3.1 - WBGene00019957
# # X 4822647 - 4834885
# lowerR = 4800000
# upperR = 4855000
# R1 <- vcf_longer %>% dplyr::filter(CHROM == "X" & POS > lowerR & POS < upperR)
# R1_hdr <- hdrs %>% 
#   dplyr::filter(chrom == "X") %>%
#   dplyr::filter(start > lowerR & start < upperR & end > lowerR & end < upperR)
# 
# R1_SV_plt <- ggplot(R1) + 
#   geom_rect(data = R1_hdr, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
#   annotate("rect", xmin = 4822647 / 1e6, xmax = 4834885 / 1e6, ymin = -Inf, ymax = Inf, fill = 'seagreen', alpha = 0.3) +
#   geom_rect(data = R1 %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = R1 %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# R1_SV_plt
# 
# 
# # C01B10.6 - WBGene00015280
# # IV 6640180 - 6645200
# lowerC = 6620000
# upperC = 6655000
# C1 <- vcf_longer %>% dplyr::filter(CHROM == "IV" & POS > lowerC & POS < upperC)
# C1_hdr <- hdrs %>% 
#   dplyr::filter(chrom == "IV") %>%
#   dplyr::filter(start > lowerC & start < upperC & end > lowerC & end < upperC)
# 
# C1_SV_plt <- ggplot(C1) + 
#   geom_rect(data = C1_hdr, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
#   annotate("rect", xmin = 6640180 / 1e6, xmax = 6645200 / 1e6, ymin = -Inf, ymax = Inf, fill = 'seagreen', alpha = 0.3) +
#   geom_rect(data = C1 %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = C1 %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# C1_SV_plt
# 
# 
# 
# # C49C3.4 - WBGene00008194
# # IV 17320795 - 17328673
# lowerC4 = 17300000
# upperC4 = 17350000
# C4 <- vcf_longer %>% dplyr::filter(CHROM == "IV" & POS > lowerC4 & POS < upperC4)
# C4_hdr <- hdrs %>% 
#   dplyr::filter(chrom == "IV") %>%
#   dplyr::filter(start > lowerC4 & start < upperC4 & end > lowerC4 & end < upperC4)
# 
# C4_SV_plt <- ggplot(C4) + 
#   geom_rect(data = C4_hdr, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
#   annotate("rect", xmin = 17320795 / 1e6, xmax = 17328673 / 1e6, ymin = -Inf, ymax = Inf, fill = 'seagreen', alpha = 0.3) +
#   geom_rect(data = C4 %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = C4 %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# C4_SV_plt
# 
# # K06G5.1 -WBGene00010605
# # X 14214973 - 14219045
# lowerK2 = 14200000
# upperK2 = 14230000
# K2 <- vcf_longer %>% dplyr::filter(CHROM == "X" & POS > lowerK2 & POS < upperK2)
# K2_hdr <- hdrs %>% 
#   dplyr::filter(chrom == "X") %>%
#   dplyr::filter(start > lowerK2 & start < upperK2 & end > lowerK2 & end < upperK2)
# 
# K2_SV_plt <- ggplot(K2) + 
#   geom_rect(data = K2_hdr, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
#   annotate("rect", xmin = 14214973 / 1e6, xmax = 14219045 / 1e6, ymin = -Inf, ymax = Inf, fill = 'seagreen', alpha = 0.3) +
#   geom_rect(data = K2 %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = K2 %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# K2_SV_plt
# 
# # gly-1 - WBGene00001626 
# # II 10899985 - 10902204 
# lowerG = 10880000
# upperG = 10920000
# G <- vcf_longer %>% dplyr::filter(CHROM == "II" & POS > lowerG & POS < upperG)
# G_hdr <- hdrs %>% 
#   dplyr::filter(chrom == "II") %>%
#   dplyr::filter(start > lowerG & start < upperG & end > lowerG & end < upperG)
# 
# G_SV_plt <- ggplot(G) + 
#   geom_rect(data = G_hdr, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf), fill = 'gray60') +
#   annotate("rect", xmin = 10899985 / 1e6, xmax = 10902204 / 1e6, ymin = -Inf, ymax = Inf, fill = 'seagreen', alpha = 0.3) +
#   geom_rect(data = G %>% dplyr::filter(SV != "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200 + SV_length)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   geom_rect(data = G %>% dplyr::filter(SV == "DEL"), aes(xmin = (POS - 200)/1e6, xmax = (POS + 200)/1e6, ymin=0.25, ymax=0.75, fill=SV)) +
#   facet_wrap(~sample) +
#   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue", "INV" = 'gold')) +
#   theme(
#     panel.background = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     panel.border = element_rect(fill = NA),
#     axis.title = element_text(size = 13),
#     legend.position = "none")
# G_SV_plt
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### OLD VISUALIZATION WITH PAFTOOLS CALLS ###
# # data <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies/elegans/77HDRs.tsv", col_names = c("CHROM","START","END","strain"))
# # 
# # all_regions <- data %>%
# #   dplyr::arrange(CHROM,START) %>%
# #   dplyr::group_split(CHROM)
# # 
# # strain_count <- data %>% dplyr::distinct(strain, .keep_all = T)
# # print(nrow(strain_count))
# # 
# # # Collapsing all HDRs
# # getRegFreq <- function(all_regions) {
# #   all_collapsed <- list()
# #   for (i in 1:length(all_regions)) {
# #     temp <- all_regions[[i]]
# #     k=1
# #     j=1
# #     while (k==1) {
# #       print(paste0("chrom:",i,"/iteration:",j))
# #       checkIntersect <- temp %>% 
# #         dplyr::arrange(CHROM,START) %>%
# #         dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
# #         dplyr::mutate(check=ifelse(is.na(check),F,check))
# #       
# #       #print(nrow(checkIntersect %>% dplyr::filter(check==T)))
# #       
# #       if(nrow(checkIntersect %>% dplyr::filter(check==T)) == 0) {
# #         print("NO MORE INTERSECTS")
# #         k=0
# #       } else {
# #         
# #         temp <- checkIntersect %>%
# #           dplyr::mutate(gid=data.table::rleid(check)) %>%
# #           dplyr::mutate(gid=ifelse((check==F| is.na(check)) & lag(check)==T,lag(gid),gid))
# #         
# #         collapse <- temp %>%
# #           dplyr::filter(check==T | (check==F & lag(check)==T)) %>%
# #           dplyr::group_by(gid) %>%
# #           dplyr::mutate(newStart=min(START)) %>%
# #           dplyr::mutate(newEnd=max(END)) %>%
# #           dplyr::ungroup() %>%
# #           dplyr::distinct(gid,.keep_all = T)  %>%
# #           dplyr::mutate(START=newStart,END=newEnd) %>%
# #           dplyr::select(-newEnd,-newStart)
# #         
# #         retain <- temp %>%
# #           dplyr::filter(check==F & lag(check)==F)
# #         
# #         temp <- rbind(collapse,retain) %>%
# #           dplyr::select(-gid,-check)
# #         
# #         j=j+1
# #       }
# #     }
# #     print(head(temp))
# #     all_collapsed[[i]] <- temp
# #   }
# #   return(all_collapsed)
# # }
# # 
# # HDR_collapse_master <- getRegFreq(all_regions)
# # 
# # all_collapsed <- ldply(HDR_collapse_master, data.frame) %>%
# #   dplyr::select(-strain) 
# # 
# # colnames(all_collapsed) <- c("chrom","start","end")
# 
# 
# # Change VCF to PAV INDEL calls
# # # vcf <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/vcf/updated_genomes/elegans.merged.1kbCOV.5kbALIGN.annotated.final.vcf") #JU2617 is not present (correctly, it is not an isotype ref so we don't have HDRs for it) 
# # # I changed PB306 to it's correct name of ECA259
# # 
# # 
# # vfilt <- vcf %>%
# #   dplyr::filter(INFO == "DEL" | INFO == "INS") %>%
# #   dplyr::filter(end - POS > 50) %>%
# #   dplyr::filter(`#CHROM` != "MtDNA") %>% 
# #   dplyr::rename(CHROM=`#CHROM`)
# # 
# # # print(colnames(vfilt))
# # 
# # pos_bed <- vfilt %>% 
# #   dplyr::select(CHROM,POS, INFO) %>%
# #   dplyr::rename(start=POS) %>%
# #   dplyr::mutate(end=start+1) %>% 
# #   dplyr::rename(chrom=CHROM)
# # 
# # 
# # 
# # ggplot() + 
# #   # geom_rect(data = all_collapsed %>% dplyr::rename(CHROM = chrom), aes(xmin = start/1e6, xmax = end/1e6, ymin = 0, ymax = 1.75, fill = "hdr")) + 
# #   geom_rect(data= vfilt %>% dplyr::filter(INFO == "DEL"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 1, ymax = 1.5, fill = "DEL")) +
# #   geom_rect(data = vfilt %>% dplyr::filter(INFO == "INS"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 0.25, ymax = 0.75, fill = "INS")) +
# #   # geom_rect(data = )
# #   facet_wrap(~CHROM, nrow=6, scales = "free_x") + 
# #   scale_fill_manual(
# #     name = "",
# #     values = c("DEL" = "red", 'INS' = 'blue', 'hdr' = 'gray70'),
# #     breaks = c("Deletions", "Insertion")
# #   ) + 
# #   theme(
# #     axis.text.y = element_blank(),
# #     axis.title.y = element_blank(),
# #     axis.ticks.y = element_blank(),
# #     axis.text.x = element_text(size = 10, color = "black"),  
# #     axis.title.x = element_text("Genome position (Mb)", size = 12, color = "black"),
# #     panel.grid = element_blank(), 
# #     plot.background = element_rect(fill = "white", color = NA),  
# #     axis.line.x = element_line(),
# #     legend.position = "none",
# #     strip.text = element_text(size = 14, color = 'black'),
# #     panel.background = element_blank(),
# #     # plot.margin = margin(0, 0, 0, 0)
# #   )
# # 
# # 
# # # isec <- valr::bed_intersect(pos_bed,all_collapsed)
# # # isec_del <- isec %>% dplyr::filter(INFO.x == "DEL")
# # # isec_ins <- isec %>% dplyr::filter(INFO.x == "INS")
# # # 
# # # 
# # # print(nrow(isec_del)) # 14747
# # # print(nrow(isec_ins)) # 17933
# # # print(nrow(isec)) # 32680
# # # 
# # # h <- nrow(isec)
# # # print(nrow(vfilt %>% dplyr::filter(INFO == "DEL"))) # 42331
# # # print(nrow(vfilt %>% dplyr::filter(INFO == "INS"))) # 59374
# # # print(nrow(vfilt)) # 101705
# # # 
# # # percent = (h/total) * 100
# # # print(percent) # 
# # 
# # 
# # vfilt_longer <- vfilt %>%
# #   dplyr::select(-end,-REF,-ALT) %>%
# #   pivot_longer(
# #     cols = -c(CHROM, POS, INFO),  # Keep chrom, pos, and info as identifiers
# #     names_to = "sample",          # New column for sample names
# #     values_to = "genotype"           # New column for sample values
# #     ) %>%
# #   dplyr::mutate(genotype=ifelse(genotype=="./.",0,1)) %>%
# #   dplyr::group_by(sample) %>%
# #   dplyr::mutate(n_alt=sum(genotype)) %>%
# #   dplyr::ungroup() %>%
# #   dplyr::arrange(desc(n_alt)) %>%
# #   dplyr::mutate(rid=rleid(n_alt)) %>%
# #   dplyr::filter(genotype > 0)
# # 
# # ### Then you will need to intersect with HDRs for a particular strain and add a column to indicate if that variant is found in a HDR for that particular strain
# # INDEL_plt <- ggplot(vfilt_longer) + geom_rect(aes(xmin=(POS-500)/1e6,xmax=(POS+500)/1e6,ymin=rid+0.7,ymax=rid-0.7,fill=INFO)) + 
# #   facet_wrap(~CHROM,nrow=1,scales = 'free_x') +
# #   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue")) +
# #   theme(
# #     panel.background = element_blank(),
# #     axis.ticks.y = element_blank(),
# #     axis.text.y = element_blank(),
# #     panel.border = element_rect(fill = NA),
# #     axis.title = element_text(size = 13),
# #     legend.position = "none") +
# #   xlab("Genome position (Mb)") +
# #   scale_y_continuous(expand = c(0,0)) + 
# #   scale_x_continuous(expand = c(0,0)) + 
# #   ylab("77 WI genomes")
# # INDEL_plt
# # 
# # 
# # # Plot only INDELs not in a HDR for each strain - code will look something like this
# # INDEL_plt_noHDR <- ggplot(vfilt_longer %>% dplyr::filter(inHDR == "NO")) + geom_rect(aes(xmin=(POS-500)/1e6,xmax=(POS+500)/1e6,ymin=rid+0.7,ymax=rid-0.7,fill=INFO)) + 
# #   facet_wrap(~CHROM,nrow=1,scales = 'free_x') +
# #   scale_fill_manual(values = c("DEL" = "red", "INS" = "blue")) +
# #   theme(
# #     panel.background = element_blank(),
# #     axis.ticks.y = element_blank(),
# #     axis.text.y = element_blank(),
# #     panel.border = element_rect(fill = NA),
# #     axis.title = element_text(size = 13),
# #     legend.position = "none") +
# #   xlab("Genome position (Mb)") +
# #   scale_y_continuous(expand = c(0,0)) + 
# #   scale_x_continuous(expand = c(0,0)) + 
# #   ylab("77 WI genomes")
# # INDEL_plt_noHDR
# # 
# # # ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/INDELs_77strains_lowerDPI.jpg", INDEL_plt, dpi = 300, width = 14, height = 10)
# # 
