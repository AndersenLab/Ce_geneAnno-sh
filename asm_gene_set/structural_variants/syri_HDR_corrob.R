library(ggplot2)
library(dplyr)
library(readr)
library(data.table)

nic_hdr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/output/CB4856_NIC_HDRs.tsv", col_names = c("chrom", "start", "end", "strain")) %>%
  dplyr::select(-strain)

syri_hdr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/SyRI/output/CB4856_SYRI_HDR.tsv", col_names = c("chrom","start","end"))

nic_dt  <- as.data.table(nic_hdr)
syri_dt <- as.data.table(syri_hdr)

setnames(nic_dt,  c("chrom","start","end"), c("chrom","start","end"))
setnames(syri_dt, c("chrom","start","end"), c("chrom","start","end"))

nic_dt[, id := .I]
syri_dt[, id := .I]

setkey(nic_dt, chrom, start, end)
setkey(syri_dt, chrom, start, end)

ov <- foverlaps(nic_dt, syri_dt, type = "any", nomatch = 0L) %>%
  dplyr::mutate(overlap = TRUE)

nic_overlap <- nic_dt %>%
  dplyr::left_join(ov %>% dplyr::select(chrom, start = i.start, end = i.end, overlap), by = c("chrom","start","end")) %>%
  dplyr::mutate(overlap = ifelse(is.na(overlap), "FALSE", overlap),
                source = "Nic", y = 2)

syri_overlap <- syri_dt %>%
  dplyr::left_join(ov %>% dplyr::select(chrom, start, end, overlap), by = c("chrom","start","end")) %>%
  dplyr::mutate(overlap = ifelse(is.na(overlap), "FALSE", overlap),
                source = "SyRI", y =1)

plt_dt <- nic_overlap %>% dplyr::bind_rows(syri_overlap)

############## Looking at overlap and unique HDR calls #################
ggplot(plt_dt) +
  # geom_hline(yintercept = 1.5, color = 'black') +
  # geom_hline(yintercept = 2, color = 'black') +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = y - 0.4, ymax = y + 0.4, fill = overlap), color = "black", size = 0.2) +
  facet_wrap(~chrom, scales = "free_x", ncol = 1) +
  scale_y_continuous(breaks = c(1,2), labels = c("SyRI","Nic")) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 16, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.x = element_text(size = 18, color = 'black')) +
  labs(fill = "Overlaps?", x = "Genomic position (Mb)")

############## Stats #################
stats_nic <- nic_overlap %>%
  dplyr::distinct() %>%
  dplyr::group_by(overlap) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup() %>% dplyr::distinct(count)

stats_syri <- syri_overlap %>%
  dplyr::distinct() %>%
  dplyr::group_by(overlap) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup() %>% dplyr::distinct(count)



# Looking at distribution of HDR call size
nic_hdr_sizes <- nic_hdr %>% dplyr::mutate(size = end - start, source = "Nic") %>% dplyr::select(size,source)

syri_hdr_sizes <- syri_hdr %>% dplyr::mutate(size = end - start, source = "SyRI") %>% dplyr::select(size,source)

plt_dist <- nic_hdr_sizes %>% dplyr::bind_rows(syri_hdr_sizes)

mean_syri <- mean(syri_hdr_sizes$size)
mean_nic <- mean(nic_hdr_sizes$size)

ggplot(plt_dist) +
  geom_boxplot(aes(y = size, x = source, fill = source)) +
  geom_jitter(aes(y = size, x = source), color = 'black', width = 0.2) +
  scale_fill_manual(values = c("Nic" = "orange", "SyRI" = "blue")) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    strip.text = element_text(size = 16, color = 'black'),
    legend.position = 'none',
    axis.text = element_text(size = 14, color = 'black'),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, color = 'black')) +
  labs(y = "Size") +
  scale_y_log10()
