library(dplyr)
library(tidyr)
library(ggplot2)
library(fuzzyjoin)
library(purrr)
library(stringr)
library(data.table)
library(cowplot)

# ======================================================================================================================================================================================== #
# Loading N0.tsv with transcripts converted to genes and plotting gene set classification #
# ======================================================================================================================================================================================== #
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0504_genes.tsv")

strainCol <- colnames(ortho_genes_dd)
strainCol_c1 <- gsub(".braker.longest.protein","",strainCol)
strainCol_c2 <- gsub(".longest.protein","",strainCol_c1)
colnames(ortho_genes_dd) <- strainCol_c2

ortho_count <- ortho_genes_dd

strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))


# ======================================================================================================================================================================================== #
# Splitting complex HOGs from 1-to-1s #
# ======================================================================================================================================================================================== #
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/116_genesOnly_strainRes.tsv", col_names = c("contig","type", "start", "end", "strand", "attributes", "strain")) 
all_genes_strain <- genes_strain %>%
  dplyr::filter(strain != "N2" | grepl("protein_coding", attributes)) %>% 
  dplyr::mutate(attributes = gsub("ID=gene:","", attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","", attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) 

ws_genes <- all_genes_strain %>% dplyr::filter(strain != "N2") %>% dplyr::select(-type)

nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/nucmer_runs/115_WI_transformed_coords_FIXED.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>%
  dplyr::select(-IDY,-LENR)

all_relations_rowid <- all_relations %>% dplyr::rename(rowid = HOG)
ortho_genes_dd_rowid <- ortho_genes_dd %>% dplyr::rename(rowid = HOG)

# Extract rowids for simple and complex HOGs
simple_rowids <- all_relations_rowid %>% 
  dplyr::filter(if_all(2:ncol(.), ~ is.na(.) | . <= 1)) %>% 
  dplyr::pull(rowid)

complex_rowids <- all_relations_rowid %>% 
  dplyr::filter(if_any(2:ncol(.), ~ . > 1)) %>% 
  dplyr::pull(rowid)

simple_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% simple_rowids)
complex_HOGS <- ortho_genes_dd_rowid %>% dplyr::filter(rowid %in% complex_rowids) ### further decompress into only complex orthogroups that have ONE N2 gene

n2_goi <- c("WBGene00019435", "WBGene00001626", "WBGene00010605", 
            "WBGene00008194", "WBGene00015280", "WBGene00019957")

n2_genes <- all_genes_strain %>%
  dplyr::filter(strain == "N2") %>%
  dplyr:: rename(chrom = contig) %>%
  dplyr::filter(str_detect(attributes, str_c(n2_goi, collapse = "|"))) %>%
  data.table::as.data.table()

# ======================================================================================================================================================================================== #
# Extracting the longest WS contig alignment for every N2 gene coordinate #
# ======================================================================================================================================================================================== #
nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, contig, WSS, WSE, L2, LENQ, strain) %>%
  data.table::as.data.table()

data.table::setkey(nucmer_ranges, chrom, start, end)
data.table::setkey(n2_genes, chrom, start, end)

n2_genes_align <- data.table::foverlaps(
  x = n2_genes,
  y = nucmer_ranges,
  type = "any" # if the start/end of any alignment is within an N2 gene
) %>%
  dplyr::rename(N2 = i.strain, n2_gene = attributes, start_aln = start, end_aln = end, start_gene = i.start, end_gene = i.end) %>%
  dplyr::select(-type, -strand)
# 10,370 N2 genes 



blah <- ggplot(n2_genes_align %>% dplyr::filter(n2_gene == "WBGene00019435")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00019435 : X")
blah


nucmer_longest <- n2_genes_align %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, contig, strain) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L2)) %>% #summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times... does this mean it is a "better" alignment?
  dplyr::ungroup() %>% ##### ADD AN ARGUMENT FOR IF TWO CONTIGS HAVE THE SAME ALINGMENT, TO CHOOSE THE LONGEST LENQ
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()


# ======================================================================================================================================================================================== #
# Removing contigs that are far away in WI coordinate system - large jumps in alignment #
# ======================================================================================================================================================================================== #
nucmer_longest_jump <- nucmer_longest %>%
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) %>% # addresses inverted alignments for tigTrim
  dplyr::arrange(St2) %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(leadDiff = lead(St2)-Et2) %>%
  dplyr::ungroup()


nucmer_longest_jumpRemoved <- nucmer_longest_jump %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(leadDiff = ifelse(is.na(leadDiff), 0, leadDiff)) %>%
  dplyr::mutate(jump = ifelse(abs(leadDiff) > 4.5E5, 1, 0)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, strain, run_id) %>%
  dplyr::mutate(gsize = n()) %>%
  dplyr::mutate(len = abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(sumlen == max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup()



# ======================================================================================================================================================================================== #
# Selecting longest alignment for contigs that align twice after filtering for jumps 
# ======================================================================================================================================================================================== #
nucmer_longest_jumpRemoved <- nucmer_longest_jumpRemoved %>%
  dplyr::mutate(row_id = dplyr::row_number())

diff <- nucmer_longest_jumpRemoved %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(dplyr::n() == 2) %>%
  dplyr::arrange(L2, .by_group = TRUE) %>%
  dplyr::mutate(diff = abs(diff(L2)[1])) %>% 
  dplyr::mutate(diff_fraction = 1 - (min(L2) / max(L2))) %>%
  dplyr::arrange(St2, .by_group = TRUE) %>%
  dplyr::mutate(jumptwo = abs(lead(St2) - Et2)) %>%
  dplyr::ungroup() ### ~55,000 rows, so 27,500 different gene-strain entries, and then 237 genes per strain, so ~1 % of n2 genes for each strain

diff_dist_plot <- ggplot(data = diff) +
  geom_histogram(aes(x = diff / 1e3), bins = 500) + 
  theme_bw() +
  xlab("Pairwise difference in alignment (kb)")
diff_dist_plot       

diff_dist_plotjumps <- ggplot(data = diff) +
  geom_histogram(aes(x = jumptwo / 1e6), bins = 200, fill = 'black') +
  geom_vline(xintercept = max(diff$jumptwo, na.rm = TRUE) / 1e6, color = 'red', size = 2) +
  # geom_vline(xintercept = 0.5, color = 'blue', size = 2) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  xlab("Difference in alignment jumps (absolute value) (Mb)")
diff_dist_plotjumps


diff_info <- diff %>%
  dplyr::mutate(n2_gene_len = (end_gene - start_gene)) %>%
  dplyr::mutate(n2_gene_middle = start_gene + (n2_gene_len / 2)) %>%
  dplyr::mutate(slope = ((Et2 - St2) / (end_aln - start_aln))) %>%
  dplyr::mutate(intercept = St2 - (slope * start_aln)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_n2_middleGene = ((slope * n2_gene_middle) + intercept)) %>%
  dplyr::arrange(St2) %>%
  dplyr::group_by(n2_gene, strain) %>% 
  dplyr::mutate(local_dup = ifelse((min(Et2) < max(St2)) & (max(St2) > min(Et2)) & ( (L1 == min(L1) & start_aln < (start_gene - 6000) & end_aln > (end_gene + 6000)) ), TRUE, FALSE)) %>%
  dplyr::mutate(WS_n2_middleGene_diff = max(WS_n2_middleGene) - min(WS_n2_middleGene)) %>% # want to include this data, because we probably don't want to keep massive jumps in duplication coordinates
  dplyr::ungroup() 


dist <- ggplot(diff_info) + 
  scale_y_continuous(expand = c(0.005, 0)) +
  scale_x_continuous(expand = c(0.005,0)) +
  geom_point(data = diff_info %>% dplyr::filter(local_dup == FALSE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_point(data = diff_info %>% dplyr::filter(local_dup == TRUE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_vline(xintercept = 100, color = "darkseagreen4", size = 2, linetype="dashed") + 
  geom_rect(xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'darkseagreen', alpha = 0.008) + 
  geom_hline(yintercept = 0.05, color = "darkseagreen4", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05, fill = 'darkseagreen', alpha = 0.008) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.title = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 14, color = 'black', face = 'bold'),
    panel.border = element_rect(fill = NA)) + 
  labs(x = "WS alignment coordinate difference at N2 gene locus (kb)", y = "Proportional difference in size of WS alignments (1 - (min(L2) / max(L2)))")
dist

# How many are local (contained in longer WS contig alignment) duplications?
local_dup_counts <- diff_info %>%
  dplyr::count(local_dup) # 175 are local, or should be kept due to syntenic context, and 232,000 do not meet the criteria

gene_jump_dist <- ggplot(data = diff_info) +
  geom_histogram(aes(x = WS_n2_middleGene_diff / 1e3), bins = 200, fill = 'gray30') + 
  scale_y_continuous(expand = c(0.001, 0)) +
  scale_x_continuous(expand = c(0.001,0)) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(r = 10, t = 10, l = 10),
    axis.line = element_line(color = "black"),           
    axis.line.y.right = element_blank(),                 
    axis.line.x.top = element_blank()) 
gene_jump_dist  

L2_diff_dist <- ggplot(data = diff_info) +
  geom_histogram(aes( y= diff_fraction), bins = 200, fill = 'gray30') + 
  scale_y_continuous(expand = c(0.001, 0), position = "right") +
  # scale_x_reverse(expand = c(0, 0.001)) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    # axis.ticks.x = element_blank(),    
    axis.text = element_text(size = 14, color = 'black'),
    panel.grid = element_blank(),
    # axis.text.x = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin( l = 20, t = 10, r = 10, b = 23),
    axis.line = element_line(color = "black"),           
    axis.line.y.right = element_blank(),                 
    axis.line.x.bottom = element_blank()) 
L2_diff_dist 


diff_filtered <- diff_info %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::mutate(
    drop_smaller = if (
      any(WS_n2_middleGene_diff > 100000 & diff_fraction > 0.05 & local_dup == FALSE)
    ) {
      (L1 != max(L1, na.rm = TRUE)) & (local_dup == FALSE)
    } else {
      FALSE  # keep both
    }
  ) %>%
  dplyr::filter(!drop_smaller) %>%
  dplyr::ungroup()

diff_filtered_twos <- diff_filtered %>%
  dplyr::group_by(n2_gene, strain) %>%
  dplyr::filter(dplyr::n() == 2) %>%
  dplyr::ungroup()

dist_filt <- ggplot(diff_filtered_twos) + 
  scale_y_continuous(expand = c(0.005, 0)) +
  scale_x_continuous(expand = c(0.005,0)) +
  geom_point(data = diff_filtered_twos %>% dplyr::filter(local_dup == FALSE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_point(data = diff_filtered_twos %>% dplyr::filter(local_dup == TRUE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_vline(xintercept = 100, color = "darkseagreen4", size = 2, linetype="dashed") + 
  geom_rect(xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'darkseagreen', alpha = 0.008) + 
  geom_hline(yintercept = 0.05, color = "darkseagreen4", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05, fill = 'darkseagreen', alpha = 0.008) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # axis.title = element_blank(),
    axis.text = element_text(size = 14, color = 'black'),
    axis.title = element_text(size = 14, color = 'black', face = 'bold'),
    panel.border = element_rect(fill = NA)) + 
  labs(x = "WS alignment coordinate difference at N2 gene locus (kb)", y = "Proportional difference in size of WS alignments (1 - (min(L2) / max(L2)))")
dist_filt

# Remove n2gene-strain pairs with two alignments from nucmer_longest_jumpRemoved, and then append filtered alignments
rows_to_remove <- diff$row_id

nucmer_longest_jumpRemoved_doublesGone <- nucmer_longest_jumpRemoved %>%
  dplyr::filter(!row_id %in% rows_to_remove)

add_to_nucmerLongest <- diff_filtered %>%
  dplyr::select(chrom,start_aln,end_aln,L1,longest_contig,WSS,WSE,L2,LENQ,strain,start_gene,end_gene,n2_gene,N2,nalign,ntig,tigsize,inv,St2,Et2,leadDiff,jump,run_id,len)

nucmer_longest_jumpRemoved_updated <- nucmer_longest_jumpRemoved_doublesGone %>%
  dplyr::bind_rows(add_to_nucmerLongest)


# ======================================================================================================================================================================================== #
# Trimming alignment(s) to ROI (N2 gene) 
# ======================================================================================================================================================================================== #
trim_spacer = 5e3 # trimming to 5kb on either side of the N2 gene
tigTrim <- nucmer_longest_jumpRemoved_updated %>%
  dplyr::arrange(n2_gene,strain,start_aln) %>% 
  dplyr::mutate(unchanged_start_aln = start_aln, unchanged_end_aln = end_aln) %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(
    scale_distortion = ((L2 - L1)/L1), # get the distortion of the WS to N2 coordinate transformation - is the slope one?
    rboundDist = max(start_aln, end_aln) - end_gene,
    Et2 = ifelse(rboundDist > trim_spacer, Et2 + (round(scale_distortion*rboundDist)) - (rboundDist - trim_spacer), Et2),
    end_aln = ifelse(rboundDist > trim_spacer,end_aln - (rboundDist - trim_spacer),end_aln),
    lboundDist = start_gene - min(start_aln, end_aln),
    St2 = ifelse(lboundDist > trim_spacer, St2 + (round(scale_distortion*lboundDist)) + (lboundDist - trim_spacer), St2),
    start_aln = ifelse(lboundDist > trim_spacer,start_aln + (lboundDist - trim_spacer),start_aln))


# ======================================================================================================================================================================================== #
# Adding WS genes to each alignment, and then collapsing by N2 gene to get a matrix of syntenic genes #
# ======================================================================================================================================================================================== #
wsg <- data.table::as.data.table(ws_genes)

clean_tigTrim <- tigTrim %>%
  dplyr::select(chrom, unchanged_start_aln, unchanged_end_aln, start_aln, end_aln, L1, start_gene, end_gene, n2_gene, WSS, WSE, inv, St2, Et2, longest_contig, nalign, L2, strain) 

filt_nucm_long <- data.table::as.data.table(clean_tigTrim)

data.table::setnames(filt_nucm_long, c("longest_contig", "St2", "Et2"), c("contig", "start", "end")) # use St2 Et2 because the gff coordinates will be reported in this way

data.table::setkey(wsg, contig, strain, start, end)
data.table::setkey(filt_nucm_long, contig, strain, start, end)

joined <- foverlaps(
  x = wsg,
  y = filt_nucm_long,
  type = "any" # if the start/end of the gene is contained within the N2 gene space or touching the boundaries of the WS alignment coordinates
) 

syntelog_matrix <- joined %>%
  dplyr::select(n2_gene,strain,attributes) %>%
  dplyr::rename(WS_gene = attributes) %>%
  tidyr::pivot_wider(
    id_cols = n2_gene,
    names_from = c(strain),
    values_from = c(WS_gene),
    values_fn = \(x) paste(unique(x), collapse = ",")
  ) %>%
  dplyr::filter(!is.na(n2_gene)) # remove the row that contains all non-syntenic predicted WS genes


######## LOOK AT GENES FOR ROBERT ######## 
# facet by strain # - also look at SV density in these gene regions + 10kb on either side? Are all of them in HDRs (shade/indicate somehow), and shade the gene in plot

# K06A9.1 - WBGene00019435
rob <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00019435")) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00019435")
rob


ab1 <- ggplot(clean_tigTrim %>% dplyr::filter(n2_gene == "WBGene00019435" & strain == "AB1")) +
  # facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00019435") +
  coord_cartesian(xlim = c(1.55,1.56), ylim = c(1.49,1.52))
ab1

ab <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00019435" & strain == "AB1")) +
  # facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 10, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00019435")
ab

# R08E3.1 - WBGene00019957
rob1 <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00019957")) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  # geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00019957")
rob1


# C01B10.6 - WBGene00015280
rob2 <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00015280")) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  # geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00015280")
rob2


# C49C3.4 - WBGene00008194
rob3 <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00008194")) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  # geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00008194")
rob3 


# K06G5.1 -WBGene00010605
rob4 <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00010605")) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  # geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00010605")
rob4

# gly-1 - WBGene00001626 
rob5 <- ggplot(joined %>% dplyr::filter(n2_gene == "WBGene00001626")) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = i.start / 1e6, ymax = i.end / 1e6), fill = "seagreen", alpha = 0.2) +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = start / 1e6, yend = end / 1e6, color = contig), linewidth = 1) +
  # geom_segment(aes(x = unchanged_start_aln / 1e6, xend = unchanged_end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = contig), linewidth = 1) +
  theme(
    legend.position = 'none',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA),
    plot.title = element_text(size = 16, color = 'black', face = 'bold', hjust = 0.5)) +
  ggtitle("WBGene00001626")
rob5

