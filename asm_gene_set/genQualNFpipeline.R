library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

merged <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/20250314_PacBio-assembly/20250314_PacBio_all_stats.txt")
single <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/merge_qual_test/WI_assembly_stats_20250321.tsv") %>%
  dplyr::filter(!(Genome_size==103831705 | Genome_size==105349000))

n50_new <- merged %>% dplyr::select(species,strain, ctg_N50, species) %>% dplyr::filter(species == "CE")
n50_old <- single %>% dplyr::select(sp,strain, N50, run, mean_read_len,yield) %>% dplyr::filter(sp == "CE")

n50 <- n50_old %>%
  dplyr::left_join(n50_new,by="strain") %>%
  dplyr::mutate(qual=ifelse(N50 <= 1e6,"bad_qual","good_qual")) %>%
  dplyr::rename(single_N50=N50,merged_N50=ctg_N50) %>%
  # dplyr::mutate(single_N50 = ifelse(is.na(single_N50), merged_N50, single_N50)) %>%
  # dplyr::mutate(run = ifelse(is.na(run), "20250305", run)) %>%
  dplyr::select(species, strain, single_N50, merged_N50, run) %>%
  dplyr::filter(!is.na(merged_N50))

lines <- n50 %>% 
  dplyr::select(strain,merged_N50,single_N50) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(g=n()) %>%
  dplyr::filter(g==2) %>%
  dplyr::mutate(x=min(merged_N50),xend=max(merged_N50),y=min(single_N50),yend=max(single_N50)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(strain,.keep_all = T) %>%
  dplyr::select(-merged_N50,-single_N50)

print(nrow(n50 %>% dplyr::distinct(strain, .keep_all = T)))

p0 <- ggplot() + 
  geom_segment(data=lines,aes(x=x/1e6,xend=xend/1e6,y=y/1e6,yend=yend/1e6),linetype='dotted',alpha=0.5) +
  geom_point(data=n50 %>% dplyr::mutate(`Sequencing Run` = as.character(run)), aes(x=merged_N50/1e6,y=single_N50/1e6, col= `Sequencing Run`), size=2.5) + 
  #geom_point(data=n50_svm,aes(x=scaf_N50/1e6,y=N50/1e6,col=class)) + 
  geom_abline(slope=1,intercept = 0) + 
  geom_hline(yintercept = 1,  linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'inside',
        legend.position.inside = c(0.2,0.8),
        axis.text = element_text(size=13),
        axis.title = element_text(size=15),
        legend.title = element_text(size = 12)) +
  xlab("Merged assembly N50 (Mb)") +
  ylab("Single assembly N50 (Mb)") 

p0

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/seq_qual_plot.png",p0, dpi = 600, height = 10, width = 15)


merged_worse <- n50 %>%
  dplyr::mutate(n50_diff = merged_N50 - single_N50)# %>%
  dplyr::select(strain,n50_diff)


  

  
  
  
  
  
  
  
data <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/FINAL_ALL_PATH_STATS_GENOMES.tsv")

data1 <- data %>%
  dplyr::filter(species == "CE") 

color_limits <- range(data1$ctg_N50 / 1e6, na.rm = TRUE)

p1 <- ggplot(data1) + 
  geom_point(aes(x=ctg_N50/1e6, y=ctg_N90, color = ctg_N50/1e6), size = 2.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")) +
  scale_color_gradient(low = "blue", high = "violet", limits = color_limits) +
  ylab("L90") + 
  xlab("N50 (Mb)") 
p1

p2 <- ggplot(data1) + 
  geom_point(aes(x=ctg_N50/1e6, y=ctg_N90, color = ctg_N50/1e6), size = 2.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")) +
  scale_color_gradient(low = "blue", high = "violet", limits = color_limits) +
  ylim(0,150) + 
  ylab("L90") + 
  xlab("N50 (Mb)")
p2

p3 <- ggplot(data1) + 
  geom_point(aes(x=fold_cov, y=ctg_N50/1e6, color = fold_cov)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")) +
  scale_color_gradient(low = "blue", high = "violet") +
  xlab("Fold Coverage") + 
  ylab("N50 (Mb)")
p3


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/n50n90_ce.png", p1, dpi = 600)
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/n50n90_ce_subset.png", p2, dpi = 600)
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/foldcovn50_ce.png", p3, dpi = 600)
