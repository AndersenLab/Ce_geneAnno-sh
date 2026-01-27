library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.table)
library(cowplot)

# ======================================================================================================================================================================================== #
# Loading in orthogroups
# ======================================================================================================================================================================================== #
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

class <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined")) %>%
  dplyr::select(Orthogroup,class)

all_class <- ortho_genes_dd %>% dplyr::bind_rows(private_OGs) %>% dplyr::left_join(class, by = "Orthogroup") 

long_class <- all_class %>%
  tidyr::pivot_longer(
    cols = -c(Orthogroup, class),
    names_to = "strain",
    values_to = "gene",
    values_drop_na = TRUE) %>%
  tidyr::separate_rows(gene, sep = ",\\s*") %>%
  dplyr::select(strain, gene, class) %>%
  dplyr::mutate(gene = sub("\\.[^.]*$", "", gene))

# pivot the table longer where each row is a gene for every strain and it's classification

# test <- all_class %>% dplyr::group_by(class) %>% dplyr::mutate(count = n()) %>% dplyr::distinct(class,count)

# Loading all genes in pangenome
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/142_140WSs_andCGC1_longestIsoGenes.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
genes_strain <- rbind(genes_strain,N2_gff)
all_genes_strain <- genes_strain %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene")


# HDRs
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
WSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/140_correct.tsv", col_names = "strain") %>% dplyr::pull()
hdrs <- hdrs %>% dplyr::filter(strain %in% WSs)

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
  dplyr::select(-strain) %>%
  data.table::as.data.table()

colnames(all_collapsed) <- c("chrom","start","end")


# ======================================================================================================================================================================================== #
# Extracting the longest WS contig alignment for every HDR #
# ======================================================================================================================================================================================== #
nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>% 
  dplyr::filter(strain != "ECA396")

nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, contig, WSS, WSE, L2, LENQ, strain) %>%
  data.table::as.data.table()

data.table::setkey(nucmer_ranges, chrom, start, end)
data.table::setkey(all_collapsed, chrom, start, end)

n2_genes_align <- data.table::foverlaps(
  x = all_collapsed,
  y = nucmer_ranges,
  type = "any" # if the start/end of any alignment is within an HDR
) %>%
  dplyr::rename(N2 = i.strain, n2_gene = attributes, start_aln = start, end_aln = end, start_gene = i.start, end_gene = i.end) %>%
  dplyr::select(-type, -strand) # 10,333 genes!



nucmer_longest <- n2_genes_align %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(strain, n2_gene) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(n2_gene, contig, strain) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
  dplyr::ungroup() %>% 
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



#### Testing how my two-variable heuristics look when I do not remove the initial jump filter
nucmer_longest_jumpRemoved <- nucmer_longest_jump
#####


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
  geom_vline(xintercept = 4.5, color = 'blue', size = 2) +
  # geom_vline(xintercept = 0.5, color = 'blue', size = 2) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  xlab("Difference in alignment jumps (absolute value) (Mb)")
diff_dist_plotjumps


# Removing alignments that are VERY small and distant from "main" alignment
# Need to find the WS coordinates that overlap with the N2 gene for each alignment, then calculate the mean (middle) WS coordinate, 
# and then if the difference is >trim_spacer (need to account for scale_distortion between WS and N2) and the smaller alignment (L1) is not at least the length of the N2 gene + 2xtrim_spacer, then drop it
# Remove rowIDs of diff dataframe from nucmer_longest_jumpRemoved, and then rbind(rows) of filtered diff column that only contains longest, most syntenic alignment
diff_info <- diff %>%
  dplyr::mutate(n2_gene_len = (end_gene - start_gene)) %>%
  dplyr::mutate(n2_gene_middle = start_gene + (n2_gene_len / 2)) %>%
  # Now taking into account Inverted alignments
  dplyr::mutate(slope = ((WSE - WSS) / (end_aln - start_aln))) %>%
  dplyr::mutate(intercept = WSS - (slope * start_aln)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_n2_middleGene = ((slope * n2_gene_middle) + intercept)) %>%
  dplyr::group_by(n2_gene, strain) %>% 
  # dplyr::mutate(local_dup = ifelse((min(WSE) < max(WSS)) & (min(WSE) < max(WSE)) & (min(WSS) < max(WSS)) & (min(WSS) < max(WSE)) & ((L1 == min(L1) & start_aln < (start_gene - 0) & end_aln > (end_gene + 0))), TRUE, FALSE)) %>% # Inversion-friendly version!
  dplyr::mutate(local_dup = ifelse(((L1 == min(L1) & start_aln < (start_gene - 0) & end_aln > (end_gene + 0))), TRUE, FALSE)) %>% # Inversion-friendly version!
  dplyr::mutate(WS_n2_middleGene_diff = max(WS_n2_middleGene) - min(WS_n2_middleGene)) %>% # want to include this data, because we probably don't want to keep massive jumps in duplication coordinates
  dplyr::ungroup() 

diff_info_noINV <- diff %>%
  dplyr::mutate(n2_gene_len = (end_gene - start_gene)) %>%
  dplyr::mutate(n2_gene_middle = start_gene + (n2_gene_len / 2)) %>%
  dplyr::mutate(slope = ((Et2 - St2) / (end_aln - start_aln))) %>%
  dplyr::mutate(intercept = St2 - (slope * start_aln)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_n2_middleGene = ((slope * n2_gene_middle) + intercept)) %>%
  dplyr::group_by(n2_gene, strain) %>% 
  dplyr::mutate(local_dup = ifelse((min(Et2) < max(St2)) & (max(St2) > min(Et2)) & ( (L1 == min(L1) & start_aln < (start_gene - 0) & end_aln > (end_gene + 0)) ), TRUE, FALSE)) %>% 
  dplyr::mutate(WS_n2_middleGene_diff = max(WS_n2_middleGene) - min(WS_n2_middleGene)) %>% # want to include this data, because we probably don't want to keep massive jumps in duplication coordinates
  dplyr::ungroup() 


dist <- ggplot(diff_info) +
  scale_y_continuous(expand = c(0.005, 0)) +
  scale_x_continuous(expand = c(0.005,0)) +
  geom_point(data = diff_info %>% dplyr::filter(local_dup == FALSE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_point(data = diff_info %>% dplyr::filter(local_dup == TRUE), aes(x = WS_n2_middleGene_diff / 1e3, y = diff_fraction, color = local_dup)) +
  geom_vline(xintercept = 100, color = "gray30", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'gray', alpha = 0.008) +
  geom_hline(yintercept = 0.05, color = "gray30", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05, fill = 'gray', alpha = 0.008) +
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
  dplyr::count(local_dup) 
# 19,459 are local and syntenic when I lower the threshold of synteny to only having to span the N2 gene, and not 6kb up- and downstream

gene_jump_dist <- ggplot(data = diff_info) +
  geom_histogram(aes(x = WS_n2_middleGene_diff / 1e3), bins = 200, fill = 'gray30') +
  scale_y_continuous(expand = c(0.001, 0), labels = scales::label_number(scale = 1e-3)) +
  scale_x_continuous(expand = c(0.001,0)) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14, color = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, color ='black'),
    # plot.margin = margin(r = 10, t = 10, l = 10),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank()) +
  ylab("Count (thousands)")
gene_jump_dist

L2_diff_dist <- ggplot(data = diff_info) +
  geom_histogram(aes( y= diff_fraction), bins = 200, fill = 'gray30') +
  scale_y_continuous(expand = c(0.001, 0), position = "right") +
  scale_x_continuous(labels = scales::label_number(scale = 1e-3)) +
  theme(
    legend.position = 'none',
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14, color = 'black'),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16, color = 'black'),
    # plot.margin = margin(l = 20, t = 10, r = 10, b = 23),
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.bottom = element_blank()) +
  xlab("Count (thousands)")
L2_diff_dist


gene_jump_dist_clean <- gene_jump_dist + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

L2_diff_dist_clean <- L2_diff_dist + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

top_row <- plot_grid(gene_jump_dist_clean, NULL, ncol = 2, rel_widths = c(0.8, 0.20))
middle_row <- plot_grid(dist, L2_diff_dist_clean, ncol = 2, rel_widths = c(0.8, 0.20)) #+ theme(plot.margin = margin(t = 20))

final_plot <- plot_grid(top_row, middle_row, nrow = 2, rel_heights = c(0.2, 0.8))

final_plot




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
  geom_vline(xintercept = 100, color = "gray30", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf, fill = 'gray', alpha = 0.008) +
  geom_hline(yintercept = 0.05, color = "gray30", size = 2, linetype="dashed") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.05, fill = 'gray', alpha = 0.008) +
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
  dplyr::select(chrom,start_aln,end_aln,L1,longest_contig,WSS,WSE,L2,LENQ,strain,start_gene,end_gene,n2_gene,N2,nalign,ntig,tigsize,inv,St2,Et2,leadDiff)

nucmer_longest_jumpRemoved_updated <- nucmer_longest_jumpRemoved_doublesGone %>%
  dplyr::bind_rows(add_to_nucmerLongest)




# ======================================================================================================================================================================================== #
# Trimming alignment(s) to HDR
# ======================================================================================================================================================================================== #
trim_spacer = 1000 # trimming to 5kb on either side of the N2 gene
tigTrim <- nucmer_longest_jumpRemoved_updated %>%
  dplyr::arrange(n2_gene,strain,start_aln) %>%
  dplyr::mutate(unchanged_start_aln = start_aln, unchanged_end_aln = end_aln) %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    scale_distortion = ((L2 - L1)/L1), # get the distortion of the WS to N2 coordinate transformation - is the slope one?
    rboundDist = max(start_aln, end_aln) - end_gene,
    WSE = ifelse(rboundDist > trim_spacer & inv == F, WSE + (round(scale_distortion*rboundDist)) - (rboundDist - trim_spacer), WSE),
    WSE = ifelse(rboundDist > trim_spacer & inv == T, WSE - (round(scale_distortion*rboundDist)) + (rboundDist - trim_spacer), WSE),
    end_aln = ifelse(rboundDist > trim_spacer,end_aln - (rboundDist - trim_spacer), end_aln),
    lboundDist = start_gene - min(start_aln, end_aln),
    WSS = ifelse(lboundDist > trim_spacer & inv == F, WSS + (round(scale_distortion*lboundDist)) + (lboundDist - trim_spacer), WSS),
    WSS = ifelse(lboundDist > trim_spacer & inv == T, WSS - (round(scale_distortion*lboundDist)) - (lboundDist - trim_spacer), WSS),
    start_aln = ifelse(lboundDist > trim_spacer, start_aln + (lboundDist - trim_spacer), start_aln)) %>%
  dplyr::ungroup()


# Confirming tigTrim is working correctly and not distorting coordinates
g10plt <- ggplot(tigTrim %>% dplyr::filter(n2_gene == "WBGene00044801")) +
  geom_rect(aes(xmin = start_gene / 1e6, xmax = end_gene / 1e6, ymin = -Inf, ymax = Inf), fill = "darkolivegreen4", alpha = 0.2) +
  facet_wrap(~strain, scales = "free", strip.position = "right") +
  geom_segment(data = nucmer_longest_jump %>% dplyr::filter(n2_gene == "WBGene00044801"), aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6), color = 'black') +
  geom_segment(aes(x = start_aln / 1e6, xend = end_aln / 1e6, y = WSS / 1e6, yend = WSE / 1e6, color = longest_contig)) +
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
  ggtitle("WBGene00044801")
g10plt


# ======================================================================================================================================================================================== #
# Adding WS genes to each alignment, and then collapsing by N2 gene to get a matrix of syntenic genes #
# ======================================================================================================================================================================================== #
wsg <- data.table::as.data.table(ws_genes)

clean_tigTrim <- tigTrim %>%
  dplyr::group_by(n2_gene,strain) %>%
  dplyr::select(chrom, unchanged_start_aln, unchanged_end_aln, start_aln, end_aln, L1, start_gene, end_gene, n2_gene, WSS, WSE, inv, longest_contig, nalign, L2, strain) %>%
  dplyr::mutate(invStartTrimmed = ifelse(inv == T, WSE, WSS), invEndTrimmed = ifelse(inv == T, WSS, WSE))

filt_nucm_long <- data.table::as.data.table(clean_tigTrim)

data.table::setnames(filt_nucm_long, c("longest_contig", "invStartTrimmed", "invEndTrimmed"), c("contig", "start", "end")) # use St2 Et2 because the gff coordinates will be reported in this way

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


