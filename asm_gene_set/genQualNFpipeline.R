library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)
library(viridis)

# merged <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/20250326_PacBio-assembly/20250326_PacBio_all_stats.txt")
# single <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/merge_qual_test/WI_assembly_stats_20250404.tsv") %>%
#   dplyr::filter(!(Genome_size==103831705 | Genome_size==105349000)) # The two individual LKC assemblies that we aren't using

merged <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/20250610_PacBio-assembly/20250610_PacBio_all_stats.txt")
single <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/merge_qual_test/WI_assembly_stats_20250609.tsv") %>%
  dplyr::filter(!(Genome_size==103831705 | Genome_size==105349000)) # The two individual LKC assemblies that we aren't using

n50_new <- merged %>% dplyr::select(species,strain, ctg_N50, species) #%>% dplyr::filter(species == "CE")
n50_old <- single %>% dplyr::select(sp,strain, N50, run, mean_read_len,yield) #%>% dplyr::filter(sp == "CE")

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
  ylab("Single assembly N50 (Mb)") +
  ggtitle("Seq Run 2025_06_10")
p0

single_better <- n50 %>%
  dplyr::filter(single_N50 > merged_N50)


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/seq_qual_plot.png",p0, dpi = 600, height = 10, width = 15)


merged_worse <- n50 %>%
  dplyr::mutate(n50_diff = merged_N50 - single_N50)# %>%
  dplyr::select(strain,n50_diff)


  

  
  
  
  
  
  
  
  
  
data <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/20251025_master_sheet.tsv")

data1 <- data %>%
  dplyr::filter(species == "CE") %>%
  dplyr::filter(strain != "ECA2888" & strain != "ECA1887" & strain != "ECA1885" & strain != "ECA2949" & strain != "ECA741") %>%
  dplyr::select(species,strain,ctg_N50,ctg_L90, fold_cov, hifi_mean_readlen) %>%
  dplyr::mutate(ctg_N50 = as.numeric(ctg_N50), ctg_L90 = as.numeric(ctg_L90), fold_cov = as.numeric(fold_cov), hifi_mean_readlen = as.numeric(hifi_mean_readlen))

n50 <- data1 %>% dplyr::filter(strain != "ECA2291")
mean_N50 <- mean(n50$ctg_N50)

color_limits <- range(data1$ctg_N50 / 1e6, na.rm = TRUE)

p1 <- ggplot(data1) + 
  geom_vline(xintercept = mean_N50 / 1e6, linetype = 'solid', color = 'black', linewidth = 2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 2) +
  geom_point(aes(x = ctg_N50 / 1e6, y = ctg_L90, 
                 color = fold_cov, size = hifi_mean_readlen)) + 
  scale_color_gradient(name = "Coverage", low = "deepskyblue", high = "red") +
  annotate("text", x = (mean_N50 / 1e6) + 7.5, y = 75, label = sprintf("Mean N50 (Mb): %.2f", mean_N50 / 1e6), hjust = 1.09, vjust = -0.3, size = 12, color = "black") +
  annotate("text", x = 6.5, y = 135, label = "1 Mb threshold", hjust = 1.09, vjust = -0.3, size = 12, color = "black") +
  scale_size_continuous(name = "Mean read length", range = c(2.5, 9.5)) + 
  theme(
    legend.position = c(0.97, 0.97),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 24),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 26, color = "black"),
    axis.title = element_text(size = 28, color = "black", face = 'bold')
  ) +
  ylab("Contig L90") + 
  xlab("Contig N50 (Mb)") 
p1

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/plots/Ce_assembly_stats_115.png", p1, height = 5, width = 9, dpi = 600)

p2 <- ggplot(data1) + 
  geom_point(aes(x=ctg_N50/1e6, y=ctg_L90, color = ctg_N50/1e6), size = 2.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")) +
  scale_color_gradient(low = "blue", high = "violet", limits = color_limits) +
  ylim(0,175) +
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


p4 <- ggplot(data1) + 
  geom_point(aes(x=mean_readlen/1e3, y=ctg_N50/1e6, color = mean_readlen)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")) +
  scale_color_gradient(low = "blue", high = "violet") +
  xlab("Mean Read Length (kb)") + 
  ylab("N50 (Mb)")
p4

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/n50n90_ce.png", p1, dpi = 600)
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/n50n90_ce_subset.png", p2, dpi = 600)
# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/foldcovn50_ce.png", p3, dpi = 600)




### Making Nx graph ###
# mamba activate seqkit
#  while read strain path; do seqkit fx2tab -l -n $path | awk -v strain=$strain '{print strain "\t" $2}'; done < all_CE_strainPath_sorted.tsv > 115_contig_lengths.tsv
contigs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/140_contigs_lengths.tsv", col_names = c("strain","contig_length")) %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(desc(contig_length), .by_group = T) %>%
  dplyr::mutate(cum_length = cumsum(contig_length)) %>%  
  dplyr::ungroup()

cgc1 <- readr::read_tsv("/vast/eande106/data/c_elegans/genomes/CGC1/CGC1_contigs_lengths.tsv", col_names = c("strain","contig_length")) %>%
  dplyr::arrange(desc(contig_length)) %>%
  dplyr::mutate(cum_length = cumsum(contig_length)) 
total_lenCGC <- cgc1 %>%
  dplyr::summarise(total_length = sum(contig_length)) %>%
  dplyr::mutate(strain = "CGC1")
cgc1_final <- cgc1 %>%
  dplyr::left_join(total_lenCGC, by = 'strain') %>%
  dplyr::mutate(cum_perc = (cum_length / total_length) * 100) %>%
  dplyr::mutate(contig_len_Mb = contig_length / 1e6)


total_length <- contigs %>%
  dplyr::group_by(strain) %>%
  dplyr::summarise(total_length = sum(contig_length)) 

contig_plotting <- contigs %>%
  dplyr::left_join(total_length, by = 'strain') %>%
  dplyr::mutate(cum_perc = (cum_length / total_length) * 100) %>%
  dplyr::mutate(contig_len_Mb = contig_length / 1e6)

left_lines <- contig_plotting %>%
  dplyr::group_by(strain) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(x_start = 0, x_end = cum_perc)

left_lines_C <- cgc1_final %>%
  dplyr::slice(1) %>%
  dplyr::mutate(x_start = 0, x_end = cum_perc)


# preview_pal <- function(pal_fun, n = 8) {
#   pal <- pal_fun(n)
#   show_col(pal)
# }
# 
# # Try:
# preview_pal(pastel_palette, 10)
# preview_pal(viridis::viridis, 10)
# preview_pal(viridis::turbo, 10)

nx <- ggplot(data = contig_plotting, aes(x = cum_perc, y = contig_len_Mb, color = strain)) + 
  geom_vline(xintercept = 50,  linetype = 'dashed', color = 'gray34', size = 1.2) + 
  geom_hline(yintercept = 1,  linetype = 'dashed', color = 'gray34', size = 1.2) +
  geom_step(direction = "vh", alpha = 0.3, size = 0.8) +
  geom_segment(data = left_lines, aes(x = x_start, xend = x_end, y = contig_len_Mb, yend = contig_len_Mb, color = strain), size = 0.8, inherit.aes = FALSE, alpha = 0.3) +
  geom_segment(data = left_lines_C, aes(x = x_start, xend = x_end, y = contig_len_Mb, yend = contig_len_Mb), color = 'gray40', size = 2, inherit.aes = FALSE) +
  geom_step(data = cgc1_final, aes(x = cum_perc, y = contig_len_Mb), direction = 'vh', color = 'gray40', size = 2) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  # geom_smooth(aes(group = 1), method = "loess", color = "black", linewidth = 1.2) +
  # scale_color_manual(values = viridis::turbo(length(unique(contig_plotting$strain)))) +
  # scale_color_manual(values = rainbow(length(unique(contig_plotting$strain)))) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    plot.margin = margin(r = 30, l = 30, t = 30, b = 30),
    panel.border = element_rect(fill = NA, color = 'black'),
    axis.title = element_text(size = 24, color = 'black', face = 'bold'),
    axis.text = element_text(size = 20, color = 'black')) + 
  labs(x = "Nx (%)", y = "Contig length (Mb)")
nx



