library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)


no_mask <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/elegans_20251008-braker_annotation/no_masked_geneCount.tsv", col_names = c("strain","count"))

soft_mask <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/elegans_softMasked_20251010-braker_annotation/soft_masked_geneCount.tsv", col_names = c("strain","masked_count")) %>%
  dplyr::filter(strain != "ECA1228" & strain != "ECA1693" & strain != "ECA730") %>%
  dplyr::left_join(no_mask, by = "strain") %>%
  tidyr::pivot_longer(c(count, masked_count), names_to = "masking", values_to = "gene_count") %>%
  dplyr::mutate(masking = recode(masking, count = "No masking", masked_count = "Soft-masked")) 

pd <- position_dodge(width = 0.8)

most <- soft_mask %>% dplyr::filter(gene_count == max(gene_count))

no  <- soft_mask %>% filter(masking == "No masking") %>% pull(gene_count)
yes <- soft_mask %>% filter(masking == "Soft-masked") %>% pull(gene_count)
meanNo  <- mean(no,  na.rm = TRUE)
meanYes <- mean(yes, na.rm = TRUE)

ggplot(soft_mask, aes(x = strain, y = gene_count, fill = masking)) +
  geom_col(position = pd, width = 0.7, color = "black") +
  # overlay lines and points exactly over the bar tops
  geom_line(aes(group = masking, color = masking), position = pd, linewidth = 0.6, show.legend = FALSE) +
  annotate("text", x = Inf, y = meanNo  + 1000, label = paste0("mean: ", round(meanNo)),
           hjust = 1.09, vjust = -0.3, size = 6, color = "black") +
  annotate("text", x = Inf, y = meanYes + 1000, label = paste0("mean: ", round(meanYes)),
           hjust = 1.09, vjust = -0.3, size = 6, color = "red") +
  annotate("text", x = 0.5, y = 25500, label = "ECA396", size = 6, hjust = -7, color = 'black') + 
  geom_hline(yintercept = meanNo, color = 'black', linewidth = 2.5) +
  geom_hline(yintercept = meanYes, color = 'red', linewidth = 2.5) +
  scale_fill_manual(values = c("No masking" = "black", "Soft-masked" = "red")) +
  scale_color_manual(values = c("No masking" = "black", "Soft-masked" = "red")) +
  labs(x = NULL, y = "Gene count", fill = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black', size = 10, angle = 75, hjust = 1),
    axis.text.y = element_text(color = 'black', size = 12),
    panel.grid = element_blank(),
    axis.title = element_text(color = 'black', size = 16)) +
  coord_cartesian(ylim = c(21000, 26000), expand = FALSE)

  