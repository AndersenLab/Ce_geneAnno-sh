library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggrepel)

preFilt_old <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/misc_old_dump/FINAL_ALL_PATH_STATS_GENOMES_GOOD_20250623.fixed.final.tsv") %>%
  dplyr::select(species, strain, genome_size, ctg_N50, ctg_L90) %>%
  dplyr::rename(Un_filtered = genome_size, ctg_N50_bef = ctg_N50, ctg_L90_bef = ctg_L90) %>%
  dplyr::mutate(strain = ifelse(strain == "PB306", "ECA259", strain)) %>%
  dplyr::filter(strain != "ECA276") #the initial sequencing was quality assembly, but wrong strain, so it will be present twice in preFilt

prefilt_julaug <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/julAug_stats.tsv") %>%
  dplyr::select(species, strain, contig_bp, ctg_N50, ctg_N90) %>%
  dplyr::rename(Un_filtered = contig_bp, ctg_N50_bef = ctg_N50, ctg_L90_bef = ctg_N90)

preFilt <- preFilt_old %>%
  dplyr::bind_rows(prefilt_julaug) %>%
  dplyr::filter(ctg_N50_bef > 1000000) %>%
  dplyr::filter(strain != "JU2617")

# prel <- preFilt %>% dplyr::filter(species == "CE") %>% dplyr::pull(strain)
# postel <- filt %>% dplyr::filter(species == "CE") %>% dplyr::pull(strain)
# setdiff(prel,postel)

  

filt <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/FINAL_allPaths_stats_20250924.tsv") %>%
  dplyr::select(strain, contig_bp, ctg_N50, ctg_L90) %>%
  dplyr::rename(Filtered = contig_bp, ctg_N50_aft = ctg_N50, ctg_L90_aft = ctg_L90) %>%
  dplyr::filter(ctg_N50_aft > 1000000) %>%
  dplyr::filter(strain != "JU2617")

all <- preFilt %>%
  dplyr::left_join(filt, by = "strain")

ce <- all %>% dplyr::filter(species == "CE") %>%
  dplyr::mutate(Un_filtered = as.numeric(Un_filtered), Filtered = as.numeric(Filtered), ctg_N50_bef = as.numeric(ctg_N50_bef), ctg_N50_aft = as.numeric(ctg_N50_aft), ctg_L90_bef = as.numeric(ctg_L90_bef),
                ctg_L90_aft = as.numeric(ctg_L90_aft))
  
cb <- all %>% dplyr::filter(species == "CB")%>%
  dplyr::mutate(Un_filtered = as.numeric(Un_filtered), Filtered = as.numeric(Filtered), ctg_N50_bef = as.numeric(ctg_N50_bef), ctg_N50_aft = as.numeric(ctg_N50_aft), ctg_L90_bef = as.numeric(ctg_L90_bef),
                ctg_L90_aft = as.numeric(ctg_L90_aft))

cn <- all %>% dplyr::filter(species == "CN") %>%
  dplyr::mutate(Un_filtered = as.numeric(Un_filtered), Filtered = as.numeric(Filtered), ctg_N50_bef = as.numeric(ctg_N50_bef), ctg_N50_aft = as.numeric(ctg_N50_aft), ctg_L90_bef = as.numeric(ctg_L90_bef),
                ctg_L90_aft = as.numeric(ctg_L90_aft))  

ct <- all %>% dplyr::filter(species == "CT") %>%
  dplyr::mutate(Un_filtered = as.numeric(Un_filtered), Filtered = as.numeric(Filtered), ctg_N50_bef = as.numeric(ctg_N50_bef), ctg_N50_aft = as.numeric(ctg_N50_aft), ctg_L90_bef = as.numeric(ctg_L90_bef),
                ctg_L90_aft = as.numeric(ctg_L90_aft))


#### C ELEGANS ####
ce_long <- ce %>%
  tidyr::pivot_longer(cols = c(Un_filtered, Filtered), names_to = "Genome", values_to = "genome_size") %>%
  dplyr::mutate(strain = factor(strain, levels = ce$strain[order(ce$Un_filtered, decreasing = FALSE)]))

gensizeCe <- ggplot(ce_long, aes(x = strain, y = genome_size / 1e6, group = strain)) +
  geom_line(linetype = "dotted") +
  geom_point(aes(color = Genome)) +
  scale_color_manual(values = c("Un_filtered" = "red", "Filtered" = "blue")) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.title = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 14, color = 'black'),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5, face = 'italic'),
    axis.text.x = element_text(size = 8, color = 'black', angle = 75, hjust = 1)) +
  xlab("Strains") +
  ylab("Genome size (Mb)") +
  ggtitle("C. elegans")
gensizeCe


ce_changed <- ce %>% dplyr::filter(ctg_N50_bef != ctg_N50_aft)

# MY16 N50 is better before filtering....
N50_ce <- ggplot(ce, aes(x = ctg_N50_bef / 1e6, y = ctg_N50_aft / 1e6)) +
  geom_point(aes(color = strain),size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_N50_bef / 1e6, yend = ctg_N50_aft / 1e6, x = ctg_N50_bef / 1e6, y = ctg_N50_bef / 1e6),linetype = "dotted", alpha = 0.6) +                    
  geom_text_repel(data = ce_changed, aes(label = strain), size = 3, max.overlaps = 20) +
   theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig N50 before (Mb)") +
  ylab("Contig N50 after (Mb)") +
  ggtitle("C. elegans")
N50_ce


ce_L90_changed <- ce %>% dplyr::filter(ctg_L90_bef < ctg_L90_aft)

# ECA2555 has an L90 of 730....
L90_ce <- ggplot(ce %>% dplyr::filter(ctg_L90_bef < 200), aes(x = ctg_L90_bef, y = ctg_L90_aft)) +
  geom_point(aes(color = strain), size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_L90_bef, yend = ctg_L90_aft, x = ctg_L90_bef, y = ctg_L90_bef, linetype = "dotted", alpha = 0.6)) +                    
  geom_text_repel(data = ce_L90_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig L90 before") +
  ylab("Contig L90 after") +
  ggtitle("C. elegans")
L90_ce





#### C BRIGGSAE ####
cb_long <- cb %>%
  tidyr::pivot_longer(cols = c(Un_filtered, Filtered), names_to = "Genome", values_to = "genome_size") %>%
  dplyr::mutate(strain = factor(strain, levels = cb$strain[order(cb$Un_filtered, decreasing = FALSE)]))

gensizeCb <- ggplot(cb_long, aes(x = strain, y = genome_size / 1e6, group = strain)) +
  geom_line(linetype = "dotted") +
  geom_point(aes(color = Genome)) +
  scale_color_manual(values = c("Un_filtered" = "red", "Filtered" = "blue")) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic'),
    axis.text.x = element_text(size = 8, color = 'black', angle = 75, hjust = 1)) +
  xlab("Strains") +
  ylab("Genome size (Mb)") +
  ggtitle("C. briggsae")
gensizeCb



cb_changed <- cb %>% dplyr::filter(ctg_N50_bef != ctg_N50_aft)

N50_cb <- ggplot(cb, aes(x = ctg_N50_bef / 1e6, y = ctg_N50_aft / 1e6)) +
  geom_point(aes(color = strain),size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_N50_bef / 1e6, yend = ctg_N50_aft / 1e6, x = ctg_N50_bef / 1e6, y = ctg_N50_bef / 1e6), linetype = "dotted", alpha = 0.6) +                    
  geom_text_repel(data = cb_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig N50 before (Mb)") +
  ylab("Contig N50 after (Mb)") +
  ggtitle("C. briggsae")
N50_cb



cb_L90_changed <- cb %>% dplyr::filter(ctg_L90_bef < ctg_L90_aft)

L90_cb <- ggplot(cb %>% dplyr::filter(ctg_L90_bef < 200), aes(x = ctg_L90_bef, y = ctg_L90_aft)) +
  geom_point(aes(color = strain), size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_L90_bef, yend = ctg_L90_aft, x = ctg_L90_bef, y = ctg_L90_bef, linetype = "dotted", alpha = 0.6)) +                    
  geom_text_repel(data = cb_L90_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig L90 before") +
  ylab("Contig L90 after") +
  ggtitle("C. briggsae")
L90_cb



#### C TROPICALIS ####
ct_long <- ct %>%
  tidyr::pivot_longer(cols = c(Un_filtered, Filtered), names_to = "Genome", values_to = "genome_size") %>%
  dplyr::mutate(strain = factor(strain, levels = ct$strain[order(ct$Un_filtered, decreasing = FALSE)]))

gensizeCt <- ggplot(ct_long, aes(x = strain, y = genome_size / 1e6, group = strain)) +
  geom_line(linetype = "dotted") +
  geom_point(aes(color = Genome)) +
  scale_color_manual(values = c("Un_filtered" = "red", "Filtered" = "blue")) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic'),
    axis.text.x = element_text(size = 8, color = 'black', angle = 75, hjust = 1)) +
  xlab("Strains") +
  ylab("Genome size (Mb)") +
  ggtitle("C. tropicalis")
gensizeCt


ct_changed <- ct %>% dplyr::filter(ctg_N50_bef != ctg_N50_aft)

N50_ct <- ggplot(ct, aes(x = ctg_N50_bef / 1e6, y = ctg_N50_aft / 1e6)) +
  geom_point(aes(color = strain),size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_N50_bef / 1e6, yend = ctg_N50_aft / 1e6, x = ctg_N50_bef / 1e6, y = ctg_N50_bef / 1e6), linetype = "dotted", alpha = 0.6) +                    
  geom_text_repel(data = ct_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig N50 before (Mb)") +
  ylab("Contig N50 after (Mb)") +
  ggtitle("C. tropicalis")
N50_ct



ct_L90_changed <- ct %>% dplyr::filter(ctg_L90_bef < ctg_L90_aft)

L90_ct <- ggplot(ct %>% dplyr::filter(ctg_L90_bef < 200), aes(x = ctg_L90_bef, y = ctg_L90_aft)) +
  geom_point(aes(color = strain), size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_L90_bef, yend = ctg_L90_aft, x = ctg_L90_bef, y = ctg_L90_bef, linetype = "dotted", alpha = 0.6)) +                    
  geom_text_repel(data = ct_L90_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig L90 before") +
  ylab("Contig L90 after") +
  ggtitle("C. tropicalis")
L90_ct




#### C NIGONI ####
cn_long <- cn %>%
  tidyr::pivot_longer(cols = c(Un_filtered, Filtered), names_to = "Genome", values_to = "genome_size") %>%
  dplyr::mutate(strain = factor(strain, levels = cn$strain[order(cn$Un_filtered, decreasing = FALSE)]))

gensizeCn <- ggplot(cn_long, aes(x = strain, y = genome_size / 1e6, group = strain)) +
  geom_line(linetype = "dotted") +
  geom_point(aes(color = Genome)) +
  scale_color_manual(values = c("Un_filtered" = "red", "Filtered" = "blue")) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text.y = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic'),
    axis.text.x = element_text(size = 8, color = 'black', angle = 75, hjust = 1)) +
  xlab("Strains") +
  ylab("Genome size (Mb)") +
  ggtitle("C. nigoni")
gensizeCn

cn_changed <- cn %>% dplyr::filter(ctg_N50_bef != ctg_N50_aft)

N50_cn <- ggplot(cn, aes(x = ctg_N50_bef / 1e6, y = ctg_N50_aft / 1e6)) +
  geom_point(aes(color = strain),size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_N50_bef / 1e6, yend = ctg_N50_aft / 1e6, x = ctg_N50_bef / 1e6, y = ctg_N50_bef / 1e6), linetype = "dotted", alpha = 0.6) +                    
  geom_text_repel(data = cn_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig N50 before (Mb)") +
  ylab("Contig N50 after (Mb)") +
  ggtitle("C. nigoni")
N50_cn
# JU1418 N50 is better before filtering.....

cn_L90_changed <- cn %>% dplyr::filter(ctg_L90_bef < ctg_L90_aft)

L90_cn <- ggplot(cn %>% dplyr::filter(ctg_L90_bef < 200), aes(x = ctg_L90_bef, y = ctg_L90_aft)) +
  geom_point(aes(color = strain), size = 2.5) +                                   
  geom_abline(slope = 1, intercept = 0) +        
  geom_segment(aes(xend = ctg_L90_bef, yend = ctg_L90_aft, x = ctg_L90_bef, y = ctg_L90_bef, linetype = "dotted", alpha = 0.6)) +                    
  geom_text_repel(data = cn_L90_changed, aes(label = strain), size = 3, max.overlaps = 20) +
  theme(
    axis.title = element_text(size = 16, color = 'black'),
    panel.background = element_blank(),
    legend.position = 'none',
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 12, color = 'black'),
    plot.title = element_text(size = 16, color = 'black', hjust = 0.5, face = 'italic')) +
  xlab("Contig L90 before") +
  ylab("Contig L90 after") +
  ggtitle("C. nigoni")
L90_cn

 
