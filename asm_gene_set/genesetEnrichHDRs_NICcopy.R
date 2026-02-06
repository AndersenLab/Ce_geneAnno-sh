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

# test <- all_class %>% dplyr::group_by(class) %>% dplyr::mutate(count = n()) %>% dplyr::distinct(class,count)

# Loading all genes in pangenome
genes_strain <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/geneAnno-nf/142_140WSs_andCGC1_longestIsoGenes.tsv", col_names = c("seqid","source", "type", "start", "end", "score", "strand", "phase", "attributes", "strain")) %>% dplyr::filter(strain != "ECA396")
N2_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3") %>% dplyr::mutate(strain="N2")
all_genes_strain <- rbind(genes_strain,N2_gff) %>%
  dplyr::mutate(attributes = gsub("ID=gene:","",attributes)) %>%
  dplyr::mutate(attributes = gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes = sub(";.*", "", attributes)) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::rename(gene = attributes)

genes_class <- all_genes_strain %>%
  dplyr::left_join(long_class, by = c("gene","strain"))


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
# Nic Moya was here

#load alignments
nucmer <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("N2S","N2E","WSS","WSE","L1","L2","IDY","LENR","LENQ","N2_chr","contig","strain")) %>% 
  dplyr::filter(strain != "ECA396") %>%
  dplyr::mutate(inv = ifelse((WSS > WSE), T, F)) %>%
  dplyr::mutate(St2 = ifelse(inv == T, WSE, WSS), Et2 = ifelse(inv == T, WSS, WSE)) # add resolution for inverted alignments - need to pull genes differently

#rename variables, arrange by WS coordinates, get lag and leading coordinates of alignments in WS space
nucmer_ranges <- nucmer %>%
  dplyr::rename(start = N2S, end = N2E, chrom = N2_chr) %>%
  dplyr::select(chrom, start, end, L1, WSS, WSE, contig, L2, LENQ, inv, strain) %>%
  dplyr::group_by(strain,chrom) %>%
  dplyr::arrange(strain,contig,WSS) %>%
  dplyr::mutate(leadS=lead(start),leadE=lead(end),leadWSS=lead(WSS),leadWSE=lead(WSE),lagS=lag(start),lagE=lag(end),lagWSS=lag(WSS),lagWSE=lag(WSE)) %>%
  dplyr::ungroup() %>%
  data.table::as.data.table()

#get HDRs
strain_hdr <- all_regions %>% 
  dplyr::rename(og_hdr_start = start, og_hdr_end = end) %>% 
  dplyr::mutate(start = ifelse(og_hdr_start >= 5000, og_hdr_start - 0, og_hdr_start), end = og_hdr_end + 0) ##################  EXTENDING HDR BOUNDARIES BY 5 KB ON BOTH SIDES TO TRY AND PULL EXTREMELY DIVERGENT REGIONS
#I use the extended coords hereafter, but you can probably set the extension to 0 and it would work anyways (do not get rid of the variables tho!)

#set keys
data.table::setkey(nucmer_ranges, strain, chrom, start, end)
data.table::setkey(strain_hdr, strain, chrom, start, end)

#find overlaps
hdr_aln <- data.table::foverlaps(
  x = strain_hdr,
  y = nucmer_ranges,
  type = "any",
  nomatch = NA) %>%
  dplyr::rename(hdr_start_extended = i.start, hdr_end_extended = i.end, N2S = start, N2E = end)  %>%
  dplyr::mutate(HDRid=paste0(strain, chrom, og_hdr_start,og_hdr_end))


#select longest contig among alingments
#select longest contig among alingments
nucmer_longest <- hdr_aln %>% # equivalent of tigFilt from haplotypePlotter.R
  dplyr::group_by(HDRid) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HDRid,contig) %>%
  dplyr::mutate(ntig = n()) %>%
  dplyr::mutate(alnsize = sum(L2)) %>% # summing the number of alignments that overlap with an N2 gene... some contigs align to a single gene many times
  dplyr::ungroup() %>% 
  dplyr::group_by(HDRid) %>%
  dplyr::filter(alnsize == max(alnsize)) %>%
  dplyr::rename(longest_contig = contig) %>%
  dplyr::filter(LENQ == max(LENQ)) %>% # to filter out alignments that are the same size, but from different contigs
  dplyr::ungroup()

clipped <- nucmer_longest %>%
  dplyr::mutate(cov_start = pmax(N2S, hdr_start_extended),
                cov_end   = pmin(N2E, hdr_end_extended)) %>%
  dplyr::filter(cov_start < cov_end)

interval_coverage <- function(starts, ends) {
  ok <- !is.na(starts) & !is.na(ends) & is.finite(starts) & is.finite(ends) & (starts < ends)
  starts <- starts[ok]
  ends   <- ends[ok]
  
  n <- length(starts)
  if (n == 0) return(0)
  
  ord <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]
  
  cur_start <- starts[1]
  cur_end   <- ends[1]
  if (n == 1) return(cur_end - cur_start)
  
  total <- 0
  for (i in 2:n) {
    if (starts[i] <= cur_end) {
      cur_end <- max(cur_end, ends[i])
    } else {
      total <- total + (cur_end - cur_start)
      cur_start <- starts[i]
      cur_end   <- ends[i]
    }
  }
  total + (cur_end - cur_start)
}

hdr_coverage <- clipped %>%
  group_by(HDRid) %>%
  summarise(
    hdr_start = first(hdr_start_extended),
    hdr_end   = first(hdr_end_extended),
    hdr_len   = hdr_end - hdr_start,
    covered_bp = interval_coverage(cov_start, cov_end),
    covered_frac = covered_bp / hdr_len,
    .groups = "drop"
  )
# 
# ggplot(hdr_coverage, aes(x = covered_frac)) + #doesn't appear to have a clear coverage split
#   geom_histogram(binwidth = 0.01) +
#   labs(x = "Fraction of HDR covered", y = "Count")

nucmer_longest <- nucmer_longest %>% dplyr::left_join(hdr_coverage %>% dplyr::select(HDRid,covered_frac),by="HDRid")    

#remove jumps by HDR alignment group
nucmer_longest_jumpRemoved <- nucmer_longest %>%
  dplyr::mutate(St2=ifelse(inv==T,WSE,WSS), Et2=ifelse(inv==T, WSS, WSE)) %>%
  dplyr::arrange(HDRid,St2) %>%
  dplyr::group_by(HDRid) %>%
  dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 7.5e4, 1 ,0)) %>% #modified to 75kb
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HDRid,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HDRid) %>%
  dplyr::filter(sumlen == max(sumlen)) %>%
  dplyr::mutate(mark_spurious = ifelse(length(unique(run_id)) > 1, TRUE, FALSE)) %>%
  dplyr::filter(!mark_spurious | dplyr::row_number() == 1) %>%
  # dplyr::select(-gsize) %>%
  dplyr::ungroup() 

# ggplot(data = nucmer_longest_jumpRemoved %>% dplyr::filter(covered_frac <0.20)) +  #doesn't appear pattern specific
#   geom_histogram(aes(x = abs(leadDiff)), binwidth = 1000, fill = 'red') +
#   theme_classic() +
#   coord_cartesian(xlim = c(0,2e5),ylim=c(0,100))

#trim boundary alignments
trim_spacer = 1e3
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
nucmer_longest_trimmed <- nucmer_longest_jumpRemoved %>%
  dplyr::rowwise() %>%
  dplyr::mutate(scale_distortion = ((L2 - L1)/L1)) %>%
  dplyr::mutate(lboundDist=hdr_start_extended-N2S) %>% #change to og start if extensions are dropped
  dplyr::mutate(rboundDist=N2E-hdr_end_extended) %>%#change to og start if extensions are dropped
  dplyr::mutate(N2E=ifelse(rboundDist>trim_spacer,(N2E-(rboundDist-trim_spacer)),N2E)) %>%
  dplyr::mutate(N2S=ifelse(lboundDist>trim_spacer,(N2S+(lboundDist-trim_spacer)),N2S)) %>%
  dplyr::mutate(WSE=ifelse(rboundDist>trim_spacer & inv==T,(WSE+(rboundDist-trim_spacer)+(rboundDist*scale_distortion)),WSE)) %>%
  dplyr::mutate(WSS=ifelse(lboundDist>trim_spacer & inv==T,(WSS-(lboundDist-trim_spacer)-(rboundDist*scale_distortion)),WSS)) %>%
  dplyr::mutate(WSE=ifelse(rboundDist>trim_spacer & inv==F,(WSE-(rboundDist-trim_spacer)-(rboundDist*scale_distortion)),WSE)) %>%
  dplyr::mutate(WSS=ifelse(lboundDist>trim_spacer & inv==F,(WSS+(lboundDist-trim_spacer)+(rboundDist*scale_distortion)),WSS)) %>%
  dplyr::ungroup()

#mark potential extensions
nucmer_mark_extend <- nucmer_longest_trimmed %>%
  dplyr::group_by(HDRid) %>%
  dplyr::mutate(WSE_extend=ifelse(inv==F & WSE == max(WSE) & N2E < hdr_end_extended, T, F),
                WSS_extend=ifelse(inv==F & WSS == min(WSS) & N2S > hdr_start_extended, T,F),
                iWSE_extend=ifelse(inv==T & WSE == min(WSE) & N2E < hdr_end_extended,T,F),
                iWSS_extend=ifelse(inv==T & WSS == max(WSS) & N2S > hdr_start_extended, T,F)) %>%
  dplyr::mutate(any_extend=ifelse(WSE_extend == T | WSS_extend == T | iWSE_extend==T | WSS_extend ==T,T,F)) %>%
  dplyr::ungroup()

#estimate extension distances
nucmer_to_extend <- nucmer_mark_extend %>% 
  dplyr::filter(any_extend==T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(extend_length_WI_lead=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE), min(leadWSS,leadWSE)-WSE,
                                             ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE), min(leadWSS,leadWSE)-WSS,NA)),
                extend_length_WI_lag=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE), WSS-max(lagWSS,lagWSE),
                                            ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE), WSE-max(lagWSS,lagWSE),NA)),
                extend_length_REF_lead=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE), leadS-N2E,
                                             ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE), N2S-leadE,NA)),
                extend_length_REF_lag=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE), N2S-lagE,
                                            ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE), lagS-N2E,NA))) %>%
  dplyr::ungroup()

#visualizes the distributions of possible extension sizes
#very similar to C.b.
#50-25kb looks good as a limit
# extendDat <- rbind(nucmer_to_extend %>% 
#                      dplyr::select(extend_length_REF_lead,extend_length_WI_lead) %>% 
#                      dplyr::rename(extend_length_WI=extend_length_WI_lead,extend_length_REF=extend_length_REF_lead),
#                    nucmer_to_extend %>% 
#                      dplyr::select(extend_length_REF_lag,extend_length_WI_lag) %>% 
#                      dplyr::rename(extend_length_WI=extend_length_WI_lag,extend_length_REF=extend_length_REF_lag)) %>%
#   dplyr::filter(!is.na(extend_length_WI) & !is.na(extend_length_REF))
# 
# sc <- ggplot(data=extendDat) + 
#   geom_point(aes(x=extend_length_REF/1e3,y=extend_length_WI/1e3),size=1) + 
#   geom_rect(xmin=0,xmax=100,ymin=-Inf,ymax=Inf,fill=NA,color="grey",linetype="dashed")+
#   geom_rect(xmin=Inf,xmax=-Inf,ymin=0,ymax=100,fill=NA,color="grey",linetype="dashed")+
#   theme_classic() + 
#   xlab("REF extension distances (kb)") + 
#   ylab("WI extension distances (kb)") +
#   scale_y_continuous(expand = c(0.01,0)) +
#   scale_x_continuous(expand = c(0.01,0))
# 
# h1 <- ggplot(data=extendDat) + 
#   geom_histogram(aes(x=extend_length_REF/1e3),binwidth = 1) + 
#   theme_classic() + 
#   xlab("") + 
#   ylab("count") + 
#   coord_cartesian(xlim=c(0,100)) +
#   scale_y_continuous(expand = c(0.01,0),
#                      labels = function(y) y / 1000,
#                      name = "count (thousand)") +
#   scale_x_continuous(expand = c(0.01,0))
# 
# h2 <- ggplot(data=extendDat) + 
#   geom_histogram(aes(y=extend_length_WI/1e3),binwidth = 1) + 
#   theme_classic() + ylab("") + 
#   coord_cartesian(ylim=c(0,100)) +
#   scale_x_continuous(labels = function(x) x / 1000,
#                      expand = c(0.01, 0),
#                      name = "count (thousand)") +
#   scale_y_continuous(expand = c(0.01,0))
# 
# 
# middle_row <- cowplot::plot_grid(
#   sc,       
#   h2, 
#   ncol = 2,
#   rel_widths = c(4, 1), 
#   align = "hv"
# )
# 
# empty_plot <- ggplot() +
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#   theme_void() 
# 
# top_row <- plot_grid(
#   empty_plot,
#   h1,
#   empty_plot,
#   ncol = 3,
#   rel_widths = c(0.09,4, 1),
#   align = "hv"
# )
# 
# final_plot <- cowplot::plot_grid(
#   top_row,      
#   middle_row,    
#   ncol = 1,
#   rel_heights = c(1, 4),  
#   align = "v"
# )
# final_plot 

# bind unmarked extensions with aligments-to-extend
nucmer_extensions <- rbind(nucmer_to_extend,nucmer_mark_extend %>% dplyr::filter(any_extend==F) %>% dplyr::mutate(extend_length_WI_lead=NA,extend_length_REF_lead=NA,extend_length_WI_lag=NA,extend_length_REF_lag=NA))

# extend alignments when conditions are met
ext_max= 5e4
nucmer_extended <- nucmer_extensions %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(N2E=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, leadS,N2E)) %>%
  dplyr::mutate(N2S=ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max,leadE ,N2S)) %>%
  dplyr::mutate(N2S=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max, lagE,N2S)) %>%
  dplyr::mutate(N2E=ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, lagS,N2E)) %>%
  dplyr::mutate(WSE=ifelse(WSE_extend==T & min(leadS,leadE) > hdr_end_extended & WSE < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, min(leadWSS,leadWSE),WSE)) %>%
  dplyr::mutate(WSS=ifelse(iWSS_extend==T & max(leadS,leadE) < hdr_start_extended & WSS < min(leadWSS,leadWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max,min(leadWSS,leadWSE) ,WSS)) %>%
  dplyr::mutate(WSS=ifelse(WSS_extend==T & max(lagS,lagE) < hdr_start_extended & WSS > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < ext_max & extend_length_REF_lag < ext_max, max(lagWSS,lagWSE),WSS)) %>%
  dplyr::mutate(WSE=ifelse(iWSE_extend==T & min(lagS,lagE) > hdr_end_extended & WSE > max(lagWSS,lagWSE) & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < ext_max & extend_length_REF_lead < ext_max, max(lagWSS,lagWSE),WSE)) %>% 
  dplyr::ungroup()

# call HDRs
WS_HDRs <- nucmer_extended %>% 
  dplyr::select(strain, longest_contig, WSS, WSE, any_extend,
                HDRid, chrom, hdr_start_extended, hdr_end_extended) %>%
  dplyr::mutate(minStart = pmin(WSS, WSE, na.rm = TRUE),
                maxEnd   = pmax(WSS, WSE, na.rm = TRUE)) %>%
  dplyr::group_by(HDRid) %>%
  dplyr::summarise(strain             = dplyr::first(strain),
                   longest_contig     = dplyr::first(longest_contig),
                   minStart           = min(minStart, na.rm = TRUE),
                   maxEnd             = max(maxEnd,   na.rm = TRUE),
                   any_extend         = dplyr::first(any_extend),
                   chrom              = dplyr::first(chrom),
                   hdr_start_extended = dplyr::first(hdr_start_extended),
                   hdr_end_extended   = dplyr::first(hdr_end_extended),
                   .groups = "drop") %>%
  dplyr::mutate(divSize=maxEnd-minStart,og_divSize=hdr_end_extended-hdr_start_extended,sizeDiff=abs(og_divSize-divSize))

##### DIAGNOSTICS ######

# diagnostic plots function
plot_hdr_workflow <- function(HDOI) {
  
  # visualize selected contig
  ctg_select_plt <- ggplot(data = nucmer_longest %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.1) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = longest_contig),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1]))
  
  # visualize pre and post jump removal
  ctg_jumprm_plt <- ggplot(data = nucmer_longest_jumpRemoved %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.3) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = longest_contig),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1], " - POST JUMP RM"))
  
  # visualize trimming after jump removal
  ctg_trim_plt <- ggplot(data = nucmer_longest_trimmed %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.3) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = longest_contig),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1], " - POST TRIM"))
  
  # visualize mark for extending
  ctg_mark_plt <- ggplot(data = nucmer_mark_extend %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.3) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = any_extend),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1], " - POST MARK"))
  
  # visualize potential extensions
  hdr_leadlag_plt <- ggplot(data = nucmer_mark_extend %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.3) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2E < hdr_end_extended),
                 aes(x = leadS / 1e6, xend = leadE / 1e6,
                     y = leadWSS / 1e6, yend = leadWSE / 1e6,
                     color = "leading"),
                 linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2S > hdr_start_extended),
                 aes(x = lagS / 1e6, xend = lagE / 1e6,
                     y = lagWSS / 1e6, yend = lagWSE / 1e6,
                     color = "lagging"),
                 linewidth = 1) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = "overlapping"),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1], " - POTENTIAL EXTENSION"))
  
  # visualize extensions
  hdr_extended_plt <- ggplot(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI[2])) +
    geom_rect(aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.3) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2E < hdr_end_extended),
                 aes(x = leadS / 1e6, xend = leadE / 1e6,
                     y = leadWSS / 1e6, yend = leadWSE / 1e6,
                     color = "leading"),
                 linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2S > hdr_start_extended),
                 aes(x = lagS / 1e6, xend = lagE / 1e6,
                     y = lagWSS / 1e6, yend = lagWSE / 1e6,
                     color = "lagging"),
                 linewidth = 1) +
    geom_segment(aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = "w/Extension"),
                 linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(HDRid == HDOI[2]),
                 aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = "overlapping"),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1], " - POST EXTENSION"))
  
  
  nS <- (WS_HDRs %>%dplyr::filter(HDRid == HDOI[2]))$minStart
  nE <- (WS_HDRs %>% dplyr::filter(HDRid == HDOI[2]))$maxEnd
  
  hdr_transformed_plt <- ggplot() +
    geom_rect(data = nucmer_extended %>% dplyr::filter(HDRid == HDOI[2]),
              mapping = aes(xmin = hdr_start_extended / 1e6, xmax = hdr_end_extended / 1e6,
                            ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.3) +
    annotate("rect",
             ymin = nS / 1e6,
             ymax = nE   / 1e6,
             xmin = -Inf,
             xmax = Inf,
             fill = "gray",
             alpha = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2E < hdr_end_extended),
                 aes(x = leadS / 1e6, xend = leadE / 1e6,
                     y = leadWSS / 1e6, yend = leadWSE / 1e6,
                     color = "leading"),
                 linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(any_extend == T) %>%
                   dplyr::filter(HDRid == HDOI[2]) %>%
                   dplyr::filter(N2S > hdr_start_extended),
                 aes(x = lagS / 1e6, xend = lagE / 1e6,
                     y = lagWSS / 1e6, yend = lagWSE / 1e6,
                     color = "lagging"),
                 linewidth = 1) +
    geom_segment(data = nucmer_extended %>%
                   dplyr::filter(HDRid == HDOI[2]),
                 aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = "w/Extension"),
                 linewidth = 1) +
    geom_segment(data = nucmer_mark_extend %>%
                   dplyr::filter(HDRid == HDOI[2]),
                 aes(x = N2S / 1e6, xend = N2E / 1e6,
                     y = WSS / 1e6, yend = WSE / 1e6,
                     color = "overlapping"),
                 linewidth = 1) +
    facet_wrap(~chrom) +
    theme(panel.border = element_rect(color = 'black', fill = NA),
          panel.background = element_blank(),
          legend.title = element_blank()) +
    labs(x = "N2 genome position (Mb)",
         y = "WS contig position (Mb)",
         title = paste0(HDOI[1], " - POST EXTENSION"))
  
  # return plot list
  list(ctg_select = ctg_select_plt,
       ctg_jumprm = ctg_jumprm_plt,
       ctg_trim   = ctg_trim_plt,
       ctg_mark   = ctg_mark_plt,
       lead_lag   = hdr_leadlag_plt,
       extended   = hdr_extended_plt,
       hdr        = hdr_transformed_plt)
}

#you can pick any test cases from this df to visualize below
View(WS_HDRs)

#set your pick here
HDOI <- c("ECA2109","ECA2109V1752300017533000") # same contig, two alignments, same size, only two so is missed by jump correction
HDOI <- c("ECA3005","ECA3005V1842400019264000") #jump correction limit is high low for this fella

HDOI <- c("ECA1493","ECA1493V25490003687000") #large, but looks legit
HDOI <- c("ECA2187","ECA2187V26820003414000") #large, but looks legit

# Lance's test cases - 3309000
HDOI <- c("ECA3088", "ECA3088II33090003360000")
HDOI <- c("ECA3088", "ECA3088X1709900017105000")
HDOI <- c("ECA3088", "ECA3088I1226600012278000")
HDOI <- c("ECA3088", "ECA3088II17030001717000")
HDOI <- c("ECA3088", "ECA3088II10400001049000")
HDOI <- c("ECA3088", "ECA3088IV1417000014189000")
HDOI <- c("ECA3088", "ECA3088X1435500014378000") # example of an extension
HDOI <- c("ECA3088", "ECA3088V1704500017120000")
HDOI <- c("ECA3088", "ECA1409V1712600017364000")
HDOI <- c("ECA723", "ECA723II1431400014341000") # largest different in N2 HDR size versus lifted-over WS HDR size
HDOI <- c("ECA2948", "ECA2948V1685600017120000") # second largest
HDOI <- c("ECA1769", "ECA1769V1515000016104000") # largest HDR among all 140 WSs
HDOI <- c("ECA1725", "ECA1725V1835700019261000") # new largest  / looks legit
HDOI <- c("ECA2151", "ECA2151V1515000016071000") # 2nd largest / checks out
HDOI <- c("ECA701", "ECA701X0288000") #wtf is this, but checks out

#get plot list for HDOI
diag_list <- plot_hdr_workflow(HDOI)

#jump correction errors visualized here
# cowplot::plot_grid(diag_list[[1]]+theme(legend.position = 'none'),diag_list[[2]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.8,1))
# 
# #trimm edge errors here
# cowplot::plot_grid(diag_list[[2]]+theme(legend.position = 'none'),diag_list[[3]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.8,1))
# 
# #mark extension errors here
# cowplot::plot_grid(diag_list[[3]]+theme(legend.position = 'none'),diag_list[[4]]+theme(axis.title.y=element_blank()),diag_list[[5]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.6,0.8,1))

#previous steps all at once - overview
cowplot::plot_grid(diag_list[[3]]+theme(legend.position = 'none'),diag_list[[4]]+theme(axis.title.y=element_blank()),diag_list[[5]],diag_list[[7]],nrow=1,align = 'h',axis = 'tb',rel_widths = c(0.6,0.8,1,1))

#find HDR clusters separated by less than 5kb
gap_clust_WS_HDRs <- WS_HDRs %>%
  dplyr::select(longest_contig,minStart,maxEnd,strain) %>%
  dplyr::rename(CHROM=longest_contig,STRAIN=strain) %>%
  dplyr::mutate(divSize=maxEnd-minStart) %>%
  dplyr::arrange(STRAIN,CHROM,minStart) %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(forGapSize=lead(minStart)-maxEnd) %>%
  dplyr::mutate(flag3g=ifelse(forGapSize<=5000,"clust","noclust")) %>%
  dplyr::mutate(dec3g=ifelse(flag3g=="clust" ,"join",
                             ifelse(flag3g=="noclust" & lag(flag3g)=="clust","join","nojoin"))) %>%
  dplyr::mutate(dec3g=ifelse(is.na(dec3g),"nojoin",dec3g)) %>%
  dplyr::ungroup()

#join the clusters
joinClust_WS_HDRs<- gap_clust_WS_HDRs %>% 
  dplyr::filter(dec3g=="join") %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::mutate(segbreak=ifelse(flag3g=="noclust",paste0(dec3g,row_number()),NA)) %>%
  tidyr::fill(segbreak,.direction = 'up') %>%
  dplyr::mutate(gid=data.table::rleid(segbreak)) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(conID=paste0(CHROM,"-",STRAIN,"-",gid)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(newStart=min(minStart),newEnd=max(maxEnd)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(conID) %>%
  dplyr::mutate(newDivSize=newEnd-newStart) %>%
  dplyr::mutate(nclust=n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(conID,.keep_all = T) %>%
  dplyr::select(-minStart,-maxEnd,-divSize) %>%
  dplyr::rename(minStart=newStart,maxEnd=newEnd,divSize=newDivSize) %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN,nclust)

#separate the unclustered regions
nojoin_WS_HDRs <- gap_clust_WS_HDRs %>%
  dplyr::group_by(STRAIN,CHROM) %>%
  dplyr::filter(!(dec3g=="join")) %>%
  dplyr::ungroup() %>%
  dplyr::select(CHROM,minStart,maxEnd,divSize,STRAIN) %>%
  dplyr::mutate(nclust=1)

#join joined and unclustered regions
#size filter
#order by divergence
all_calls_WS_HDRs<- rbind(joinClust_WS_HDRs,nojoin_WS_HDRs) %>%
  dplyr::filter(divSize/1e3 >= 5) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::rename(contig=CHROM,strain=STRAIN,strain_order=ystrain) %>%
  dplyr::select(-nclust,-ncalls,-sorter,-rleID) %>%
  dplyr::arrange(strain_order,contig,minStart)

View(all_calls_WS_HDRs) #strains are orderered from least-to-most divergent based on HDR count
