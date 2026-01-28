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
  dplyr::filter(type == "gene") %>%
  dplyr::rename(gene = attributes)

genes_class <- all_genes_strain %>%
  dplyr::left_join(long_class, by = "gene","strain")


# HDRs
hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))
WSs <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/assembly-nf/all_assemblies_sheet/140_correct.tsv", col_names = "strain") %>% dplyr::pull()

# COLLAPSING HDRS AMONG 140 WSs
all_regions <- hdrs %>%
  dplyr::filter(strain %in% WSs) %>%
  dplyr::arrange(chrom,start) %>%
  data.table::as.data.table()


# ======================================================================================================================================================================================== #
# Extracting the longest WS contig alignment for every HDR #
# ======================================================================================================================================================================================== #
nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>% 
  dplyr::filter(strain != "ECA396") %>%
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) # add resolution for inverted alignments - need to pull genes differently

nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, WSS, WSE, contig, L2, LENQ, inv, strain) %>%
  data.table::as.data.table()

# Finding overlap of WS alignments with HDRs - for an INDIVIDUAL strain:
### Once finalalized, create a function with the analysis workflow and iterate over wild strains: for (i in WSs)
SOI <- "ECA3088"

nucmer_ranges <- nucmer_ranges %>% dplyr::filter(strain == SOI)
strain_hdr <- all_regions %>% dplyr::filter(strain == SOI) %>% 
  dplyr::rename(og_hdr_start = start, og_hdr_end = end) %>% 
  dplyr::mutate(start = ifelse(og_hdr_start >= 5000, og_hdr_start - 5000, og_hdr_start), end = og_hdr_end + 5000) ##################  EXTENDING HDR BOUNDARIES BY 5 KB ON BOTH SIDES TO TRY AND PULL EXTREMELY DIVERGENT REGIONS
  
data.table::setkey(nucmer_ranges, chrom, start, end)
data.table::setkey(strain_hdr, chrom, start, end)

hdr_aln <- data.table::foverlaps(
  x = strain_hdr,
  y = nucmer_ranges,
  type = "any" # if the start/end of any alignment is within an HDR
  ) %>%
  dplyr::rename(hdr_start_extended = i.start, hdr_end_extended = i.end, N2S = start, N2E = end) %>%
  dplyr::select(-i.strain) # 10,333 genes!

num_hdrs <- nrow(strain_hdr)
num_hdr_aln <- nrow(hdr_aln %>% dplyr::distinct(hdr_start_extended))
strain_id <- hdr_aln %>% dplyr::distinct(strain) %>% dplyr::pull()
print(paste0(num_hdr_aln, " / ", num_hdrs, " HDRs have overlapping ", strain_id, " alignments"))

# Test plot
ggplot(data = hdr_aln) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))

ggplot(data = hdr_aln %>% dplyr::filter(chrom == "IV" & og_hdr_start > 12000000 & og_hdr_end < 14000000)) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


nucmer_longest <- hdr_aln %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(og_hdr_start) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(og_hdr_start, contig) %>%
  dplyr::mutate(ntig = n()) %>%
  dplyr::mutate(tigsize = sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
  dplyr::ungroup() %>% 
  dplyr::group_by(og_hdr_start) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::group_by(og_hdr_start) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()

ggplot(data = nucmer_longest) +
  geom_rect(aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.1) +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


# Calculating the WS coordinates that directly overlap with the N2 HDR boundaries
nucmer_slope <- nucmer_longest %>%
  # dplyr::group_by(og_hdr_start) %>%
  # dplyr::mutate(longest_L2 = max(L2)) %>% # filtering for the longest alignment for each HDR to calculate the slope
  # dplyr::ungroup() %>%
  dplyr::mutate(slope = ((WSE - WSS) / (N2E - N2S))) %>%
  dplyr::mutate(intercept = WSS - (slope * N2S)) %>% # to find the y-intercept using point-slope form (b = y1 - mx1)
  dplyr::mutate(WS_hdr_start = ((slope * og_hdr_start) + intercept)) %>%
  dplyr::mutate(WS_hdr_end = ((slope * og_hdr_end) + intercept)) %>%
  dplyr::group_by(og_hdr_start) %>%
  dplyr::mutate(WS_hdr_start_min = min(WS_hdr_start)) %>% # conditional based on if alignment is INV and min WSS! - what if there is only an alignment for one HDR boundary???
  dplyr::mutate(WS_hdr_end_max = max(WS_hdr_end)) %>% # conditional based on if alignment is INV and max WSE! - what if there is only an alignment for one HDR boundary???
  dplyr::ungroup() 

test <- nucmer_slope %>% dplyr::filter(chrom == "IV" & og_hdr_start > 12600000 & og_hdr_end < 13000000)
ggplot(test) + 
  geom_rect(data = test %>% dplyr::distinct(og_hdr_start, .keep_all = T), aes(xmin = og_hdr_start / 1e6, xmax = og_hdr_end / 1e6, ymin = -Inf, ymax = Inf), fill = "#DB6333", alpha = 0.3) +
  geom_rect(data = test %>% dplyr::distinct(WS_hdr_start_min, .keep_all = T), aes(xmin = -Inf, xmax = Inf, ymin = WS_hdr_start_min / 1e6, ymax = WS_hdr_end_max / 1e6), fill = 'blue', alpha = 0.3) +
  geom_segment(aes(x = og_hdr_start / 1e6, xend = og_hdr_start / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = og_hdr_end / 1e6, xend = og_hdr_end / 1e6, y = -Inf, yend = Inf), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_start_min / 1e6, yend = WS_hdr_start_min / 1e6), color = 'black') +
  geom_segment(aes(x = -Inf, xend = Inf, y = WS_hdr_end_max / 1e6, yend = WS_hdr_end_max / 1e6), color = 'black') +
  geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6, y = WSS / 1e6, yend = WSE / 1e6,  color = longest_contig), size = 1) +
  facet_wrap(~chrom) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    panel.background = element_blank()
  ) +
  labs(x = "N2 genome position (Mb)", y = "WS contig position (Mb)", title = paste0(hdr_aln$strain))


ws_n2_hdr_region <- nucmer_slope %>% 
  dplyr::group_by(og_hdr_start) %>%
  dplyr::mutate(hdr_size = og_hdr_end - og_hdr_start, ws_hdr_size = WS_hdr_end_max - WS_hdr_start_min) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(og_hdr_start, .keep_all = T)

# Shows that HDRs are very structurally different than N2 and not just high SNV density
ggplot(data = ws_n2_hdr_region) +
  geom_point(aes(x = hdr_size / 1e6, y = abs(ws_hdr_size / 1e6)), color = 'firebrick') +
  geom_line(data = data.frame(x = c(0, 1.5)), aes(x = x, y = x), linetype = "dashed") +  
  scale_x_log10(limits = c(0.009, 20)) +
  scale_y_log10(limits = c(0.009, 20)) +
  theme_bw() +
  labs(x = "N2 HDR size (Mb) log scale", y = paste0(strain_id, " HDR size (Mb) log scale"))




# ======================================================================================================================================================================================== #
# Removing contigs that are far away in WI coordinate system - large jumps in alignment #
# ======================================================================================================================================================================================== #
nucmer_longest_jump <- nucmer_longest %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) %>% # addresses inverted alignments for tigTrim
  dplyr::arrange(St2) %>%
  dplyr::group_by(hdr_start) %>%
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


