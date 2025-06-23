library(plyr)
library(readr)
library(tidyr)
library(org.Ce.eg.db)
library(GO.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(clusterProfiler) ## BiocManager::install("clusterProfiler") # need this and the next package???
library(enrichplot)
library(data.table)
library(cowplot)


ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/tropicalis/processed_data/N0_genes.tsv") %>%
  dplyr::rename(NIC58 = NIC58.ONLYPC.longest.protein, N2 = c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein)

strainCol <- colnames(ortho_genes_dd)
colnames(ortho_genes_dd) <- strainCol

ortho_count <- ortho_genes_dd

strainCol_c2_u <- strainCol[!strainCol %in% c("OG", "HOG", "Gene Tree Parent Clade")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))



# ======================================================================================================================================================================================== #
# GO - extract any-to-ones between NIC58 and N2 and lift over GO terms #
# ======================================================================================================================================================================================== #
# Extracting NIC58 - N2 any-to-ones
qn <- all_relations %>%
  dplyr::filter(NIC58_count >= 1) %>%
  dplyr::filter(N2_count == 1) 

oneoneqn <- ortho_genes_dd %>%
  dplyr::filter(HOG %in% qn$HOG) %>%
  dplyr::select(NIC58, N2)


## Extracting GO terms for N2 genes
gene_ids <- keys(org.Ce.eg.db, keytype = "WORMBASE")
tran_gene <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/input/WBGene_gene_sets/longest_iso_tranName_WBGeneID.tsv", col_names = c("transcript","gene"))

n2_wbGenes <- tran_gene %>%
  dplyr::select(gene) %>%
  dplyr::distinct(gene)

# Documentation: https://www.bioconductor.org/packages/release/data/annotation/manuals/org.Ce.eg.db/man/org.Ce.eg.db.pdf
gene2go <- AnnotationDbi::select(
  org.Ce.eg.db,
  keys = gene_ids,
  columns = c("GO", "ONTOLOGY"),
  keytype = "WORMBASE"
)
head(gene2go)


n2_gene_database <- gene2go %>%
  dplyr::select(WORMBASE, GO, ONTOLOGY) %>%
  dplyr::mutate(n2_gene = ifelse(WORMBASE %in% n2_wbGenes$gene, TRUE, FALSE)) %>%
  dplyr::mutate(n2_gene = ifelse(n2_gene == TRUE, WORMBASE, n2_gene)) 

merged  <- n2_gene_database %>%
  dplyr::select(n2_gene,GO) %>%
  dplyr::filter(n2_gene != FALSE) %>% 
  dplyr::filter(!is.na(GO)) %>%
  dplyr::rename(gene = n2_gene) # 14,033 N2 genes


# Lift over GO terms from N2 to NIC58
GO_trop <- oneoneqn %>%
  tidyr::separate_rows(NIC58, sep = ",\\s*") %>% 
  dplyr::left_join(merged, by = c("N2" = "gene"), relationship = "many-to-many") %>%
  dplyr::filter(!is.na(GO)) # 10,970 NIC58 genes


# ======================================================================================================================================================================================== #
# Lift over GO IDs from RBBH output #
# ======================================================================================================================================================================================== #
rbbh <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/tropicalis/RBBH/N2_NIC58_RBBH.tsv", col_names = c("N2","NIC58_rbbh"))
sumdf <- ortho_genes_dd %>%
  dplyr::select(HOG,NIC58, N2) %>%
  dplyr::left_join(GO_trop, by = "N2") %>%
  dplyr::select(-NIC58.y) %>%
  dplyr::rename(NIC58 = NIC58.x) %>%
  tidyr::separate_rows(N2, sep = ",") %>%
  dplyr::left_join(rbbh, by = "N2") 

rbbh_GO <- sumdf %>%
  dplyr::filter(is.na(GO) & !is.na(NIC58_rbbh)) %>%
  dplyr::select(HOG,NIC58,NIC58_rbbh,GO,N2) %>%
  dplyr::left_join(merged, by = c("N2" = "gene")) %>%
  dplyr::rename(GO_rbbh = GO.y)

final_GO <- sumdf %>%
  dplyr::left_join(rbbh_GO, by = "N2") %>%
  dplyr::select(NIC58.x, N2, GO, NIC58_rbbh.x, GO_rbbh) %>%
  dplyr::rename(NIC58 = NIC58.x, NIC58_rbbh = NIC58_rbbh.x) %>%
  dplyr::filter(!is.na(NIC58) & !is.na(N2)) %>%
  tidyr::separate_rows(NIC58, sep = ",\\s") %>%
  dplyr::mutate(GO_final = ifelse(!is.na(GO), GO, GO_rbbh)) %>%
  dplyr::select(NIC58, N2, GO_final) %>%
  dplyr::filter(!is.na(GO_final)) %>%
  dplyr::select(NIC58, GO_final) %>%
  dplyr::rename(GO = GO_final) # 11,635 NIC58 genes


# ======================================================================================================================================================================================== #
# GO enrichment analysis #
# ======================================================================================================================================================================================== #
classifyingArms <- function(HDRfile, species_armDomain)
{
  HDRfile <- HDRfile %>%
    dplyr::rename(C=CHROM, S=minStart, E=maxEnd) %>%
    dplyr::select(C, S, E)
  
  df_list <- list() # set empty list to append lists of HDRs found within
  # chromosomal arms
  
  for (i in 1:nrow(species_armDomain)) {
    domChrom = as.character(species_armDomain[i,1]) # iterates by chromosome
    domStart = as.numeric(species_armDomain[i,2]) #  1e3
    domEnd = as.numeric(species_armDomain[i,3]) #* 1e3
    
    if (i %% 2 == 0){ # if i is even (the iteration number)
      D="R"
    } else { # if i is odd
      D="L"
    }
    HDRfiltered <- HDRfile %>%
      dplyr::filter(C == domChrom) %>%
      dplyr::filter((E >= domStart & E <= domEnd) | (S >= domStart & S <= domEnd)) %>%
      dplyr::mutate(domain = paste0(D, "_", "arm", "_", domChrom))
    
    df_list[[i]] <- HDRfiltered
  }
  
  species_HDRs_wDomains <- ldply(df_list, data.frame) # creating dataframe from
  # lists of HDRs found within chromosome arms
  
  nonHDR_arms <- species_HDRs_wDomains %>%
    dplyr::arrange(C, S) %>%
    dplyr::group_by(domain) %>%
    dplyr::mutate(newEnd = (lead(S)-1)) %>%
    dplyr::mutate(newStart = (E+1)) %>%
    dplyr::select(C,domain,newStart,newEnd) %>%
    dplyr::ungroup() # isolating regions found between HDRs in chrom arms
  
  start_tips <- species_HDRs_wDomains %>%
    dplyr::group_by(domain) %>%
    dplyr::filter(S == min(S) & (!is.na(E))) %>%
    dplyr::mutate(newEnd = S-1)
  
  start_tips$newStart <- species_armDomain$left #* 1e3 # adding chrom arm tip regions starts
  
  start_tips_final <- start_tips %>%
    dplyr::mutate(dropMark=ifelse(S < newStart,'D','ND')) %>% #if the start of the HD region extends beyond the arm boundary, we flag the row
    dplyr::filter(dropMark=="ND") %>% #we filter out the flagged rows
    dplyr::select(C,domain,newStart,newEnd)
  
  end_tips <- nonHDR_arms %>%
    dplyr::filter(is.na(newEnd))
  
  
  
  
  end_tips$newEnd <- species_armDomain$right #* 1e3 # adding chrom arm tip regions ends
  
  #similar error is found in the end tips, commented edits below
  end_tips_final <- end_tips  %>%
    dplyr::mutate(dropMark=ifelse(newStart>newEnd,"D","ND"))  %>% # we flag end tip rows that extend past the domain end boundary
    dplyr::filter(dropMark=='ND') %>% #and we drop those flagged rows
    dplyr::select(-dropMark)
  
  nonHDR_arms_temp <- nonHDR_arms %>%
    dplyr::filter(!is.na(newEnd))
  
  nonHDR_arms_final <- rbind(nonHDR_arms_temp, end_tips_final, start_tips_final) %>% # finalizing dataframe with tips and ends added
    dplyr::arrange(C,domain,newStart)
  
  colnames(nonHDR_arms_final) <- c("CHROM","domain","start","end")
  return(nonHDR_arms_final)
}

### Collapsing HDRs found among all supplies strains ###
getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      #print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>%
        dplyr::arrange(CHROM,minStart) %>%
        dplyr::mutate(check=ifelse(lead(minStart) <= maxEnd,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
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
          dplyr::mutate(newStart=min(minStart)) %>%
          dplyr::mutate(newEnd=max(maxEnd)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(minStart=newStart,maxEnd=newEnd) %>%
          dplyr::select(-newEnd,-newStart)
        
        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print("WRITING TO LIST")
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}


getGenesFromBed <- function(gff,regions) {
  regGenes <- list()
  for (i in 1:nrow(regions)) {
    #print(i)
    bs=regions[i,]$start
    be=regions[i,]$end
    bc=regions[i,]$CHROM
    
    tmp <- gff %>%
      dplyr::filter(CHROM==bc) %>%
      dplyr::filter((start >= bs & start <= be) &
                      (end >= bs & end <= be))
    regGenes[[i]] <- tmp
  }
  return(regGenes)
}


# nHDR_armDomain_BRIG <- readr::read_csv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/briggsae/multi_reference/2024_release/QX1410/chromosome_domain_Cbriggsae.csv") %>%
#   dplyr::rename(Chromosome=chrom) %>%
#   dplyr::select(Chromosome,left,right) %>%
#   dplyr::mutate(left=left*1e3,right=right*1e3)

nHDR_armDomain <- readr::read_tsv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/tropicalis/arm_domain_ranges.tsv") %>%
  dplyr::rename(left=domStart,right=domEnd) %>%
  dplyr::select(Chromosome,left,right)

hdr_regions <- readr::read_tsv("/vast/eande106/data/c_tropicalis/WI/divergent_regions/20250627/20250627_c_tropicalis_divergent_regions_strain.bed", col_names = c("CHROM","minStart","maxEnd","STRAIN")) 

HDreg_all <- ldply(getRegFreq(hdr_regions %>% dplyr::select(CHROM,minStart,maxEnd,STRAIN) %>% dplyr::arrange(CHROM,minStart) %>%dplyr::group_split(CHROM)), data.frame) %>% dplyr::select(-STRAIN) %>% dplyr::rename(start=minStart,end=maxEnd)

centers <- nHDR_armDomain %>%
  dplyr::group_by(Chromosome) %>%
  dplyr::mutate(center_end=lead(left)) %>%
  dplyr::filter(!is.na(center_end)) %>%
  dplyr::select(-left) %>%
  dplyr::rename(center_start=right)

reglist <- list()
reglist_center <- list()
for (i in 1:nrow(centers)) {
  cstart <- centers[i,]$center_start
  cend <- centers[i,]$center_end
  cchrom <- centers[i,]$Chromosome
  temp <- HDreg_all %>%
    dplyr::filter(CHROM==cchrom) %>%
    dplyr::filter(start < cstart | end > cend)
  temp2 <- HDreg_all %>%
    dplyr::filter(CHROM==cchrom) %>%
    dplyr::filter(!start < cstart & !end > cend)
  reglist[[i]] <- temp
  reglist_center[[i]] <- temp2
}

# Read in GFF for extracting NIC58 genes in nHDRs and HDRs
gffCt <-  ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/tropicalis/prep/NIC58.update.April2025.noWBGeneID.csq.ONLYPC.longest.gff3")

gffFinal <- gffCt %>%
  dplyr::filter(type=="gene") %>%
  dplyr::filter(grepl("biotype=protein_coding", attributes)) %>%
  tidyr::separate(attributes, into=c("ID","rest"), sep = "ID=gene:") %>%
  tidyr::separate(rest, into=c("keep","nope"), sep = ";biotype=") %>%
  dplyr::select(-ID,-source,-type,-score,-strand,-phase,-nope) %>%
  dplyr::rename(CHROM = seqid, NIC58 = keep) %>%
  dplyr::distinct(NIC58, .keep_all = T)


HDreg_arm <- ldply(reglist,data.frame)
HDreg_center <- ldply(reglist_center,data.frame)
HDgene <- getGenesFromBed(gffFinal,HDreg_arm)
HD_gene_vector <- HDgene[[]]
HD_NIC58_genes <- (do.call(rbind, HDgene))
HD_gene_vector <- unique(HD_NIC58_genes$NIC58) 

nHDreg <- classifyingArms(HDreg_arm %>% dplyr::rename(minStart=start,maxEnd=end), nHDR_armDomain) %>% dplyr::select(-domain)
nHDgene <- getGenesFromBed(gffFinal,nHDreg)
nHD_NIC58_genes <- do.call(rbind, nHDgene)
nHD_gene_vector <- unique(nHD_NIC58_genes$NIC58) 


posplot <- ggplot() + 
  geom_rect(data=centers %>% dplyr::rename(CHROM=Chromosome), aes(xmin=center_start,xmax=center_end,ymin=-1,ymax=1,fill="nHDR_center")) +
  geom_rect(data=rbind(HDreg_arm %>% dplyr::mutate(class="HDR_arm"),nHDreg %>% dplyr::mutate(class="nHDR_arm"), HDreg_center %>% dplyr::mutate(class="HDR_center")), aes(xmin=start,xmax=end,ymin=-1,ymax=1,fill=class)) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values=c("nHDR_center"="lightblue","nHDR_arm"="deepskyblue3","HDR_center"="pink","HDR_arm"="firebrick3"))+
  xlab("Physical postion (Mb)")
posplot

posplot2 <- ggplot() + 
  geom_rect(data=HD_NIC58_genes, aes(xmin=start,xmax=end,ymin=-1,ymax=1)) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  xlab("Physical postion (Mb)")
posplot2



NIC58_merged_bgd_genes <- unique(final_GO$NIC58) #  11,635 genes
GO_terms_merged <- AnnotationDbi::select(GO.db, 
                                         keys=unique(final_GO$GO),   
                                         columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                         keytype="GOID") %>% 
  dplyr::rename(TERM = GOID, TERM_NAME = TERM) # 6,772

merged_ont <- final_GO %>%
  dplyr::left_join(GO_terms_merged, by = c("GO" = "TERM")) 

# HDR genes in arms against all NIC58 background set
enGO_HDR_merged <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,NIC58),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = NIC58_merged_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged)

dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")

enGO_HDR_dt <- as.data.table(enGO_HDR_merged@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())


insert_linebreaks <- function(s, width = 40) {
  sapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

GO_list_BP <- enGO_HDR_dt %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>% 
  # dplyr::slice_head(n=30) %>%
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number())

# GO_list_BP <- enGO_HDR_dt %>%
#   dplyr::filter(p.adjust < 0.05) %>%
#   dplyr::arrange(desc(p.adjust)) %>%  
#   dplyr::mutate(Description = factor(Description, levels = unique(Description)), plotpoint = dplyr::row_number())

# Dimensions for figure #
# plot_GO_BP <- ggplot(GO_list_BP) +
#   geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
#   geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
#   scale_y_continuous(breaks = GO_list_BP$plotpoint, labels = GO_list_BP$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
#   scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(GO_list_BP$Count, na.rm = TRUE)), round((max(GO_list_BP$Count, na.rm = TRUE) + min(GO_list_BP$Count, na.rm = TRUE) ) / 2), round(max(GO_list_BP$Count, na.rm = TRUE)))) +
#   scale_size_continuous(range = c(0.1, 3), name = "Fold enrichment", breaks = pretty(GO_list_BP$EnrichRatio, n = 3)) +
#   theme(axis.text.x = element_text(size=6, color='black'),
#         axis.text.y = element_text(size=6, color='black'),
#         axis.title = element_text(size=7, color='black', face = 'bold'),
#         # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         plot.title = element_blank(),
#         legend.title = element_text(size=5, color='black', hjust = 1),
#         legend.text = element_text(size=4, color='black', hjust = 1),
#         legend.position = "inside",
#         legend.position.inside = c(0.8, 0.195),
#         legend.direction = "horizontal", legend.box = "vertical",
#         legend.spacing.y = unit(0.0001, 'cm'),
#         legend.key.height = unit(0.01, "cm"),
#         legend.key.width = unit(0.5, "cm"),
#         legend.box.just = "right",
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(fill = NA),
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
#   # legend.spacing.x = unit(0.02, "cm"),
#   # legend.key.size = unit(0.3, 'cm')) +
#   guides(
#     fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
#     size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
#   labs(title = "Enriched GO:BP terms for genes in HDRs", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
# plot_GO_BP


plot_GO_BP <- ggplot(GO_list_BP) +
  geom_vline(xintercept = -log10(0.05), color='red', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_BP$plotpoint, labels = GO_list_BP$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "skyblue", high = "#0719BC", breaks = c(round(min(GO_list_BP$Count, na.rm = TRUE)), round((max(GO_list_BP$Count, na.rm = TRUE) + min(GO_list_BP$Count, na.rm = TRUE) ) / 2), round(max(GO_list_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(GO_list_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=12, color='black'),
        axis.text.y = element_text(size=12, color='black'),
        axis.title = element_text(size=12, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        # plot.title = element_blank(),
        legend.title = element_text(size=10, color='black', hjust = 1),
        legend.text = element_text(size=8, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.195),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.01, 'cm'),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 1.5),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDR arms", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP

# GOBP_tableSave <- enGO_HDR_dt %>%
#   dplyr::filter(p.adjust < 0.05) %>%
#   dplyr::arrange(desc(p.adjust)) %>% 
#   dplyr::select(-plotpoint) %>%
#   dplyr::mutate(geneID = gsub("/", ", ", geneID))


# Need to revert-back any simplifed Descritpions
# write.table(GOBP_tableSave, file = "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOBP_HDRgenes.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,NIC58),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = NIC58_merged_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")

enGO_HDR_merged_MF <- as.data.table(enGO_HDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

insert_linebreaks <- function(s, width = 40) {
  sapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

GO_list_MF <- enGO_HDR_merged_MF %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>% 
  # dplyr::slice_head(n=30) %>%
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number()) 

# GO_list_MF <- GO_list_MF %>%
#   dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen","oxidoreductase activity (1)", Description), 
#                 Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (2)", Description),
#                 Description = gsub("transmitter-gated ion channel activity involved in regulation of postsynaptic membrane potential", "transmitter-gated ion channel (1)", Description))


# oxidoreductase activity (1) = oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen
# oxidoreductase activity (2) = oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen
# transmitter-gated ion channel (1) = transmitter-gated ion channel activity involved in regulation of postsynaptic membrane potential


# ### Dimensions for figure ###
# plot_GO_MF <- ggplot(GO_list_MF) +
#   geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
#   geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
#   scale_y_continuous(breaks = GO_list_MF$plotpoint, labels = GO_list_MF$Description, name = "", expand = c(0.02,0.02)) +
#   scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(GO_list_MF$Count, na.rm = TRUE)), round((max(GO_list_MF$Count, na.rm = TRUE) + min(GO_list_MF$Count, na.rm = TRUE) ) / 2), round(max(GO_list_MF$Count, na.rm = TRUE)))) +
#   scale_size_continuous(range = c(0.1, 3), name = "Fold enrichment", breaks = pretty(GO_list_MF$EnrichRatio, n = 3)) +
#   theme(axis.text.x = element_text(size=6, color='black'), 
#         axis.text.y = element_text(size=6, color='black'),
#         axis.title = element_text(size=7, color='black', face = 'bold'), 
#         # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         plot.title = element_blank(),
#         legend.title = element_text(size=5, color='black', hjust = 1), 
#         legend.text = element_text(size=4, color='black', hjust = 1), 
#         legend.position = "inside", 
#         legend.position.inside = c(0.8, 0.195),
#         legend.direction = "horizontal", legend.box = "vertical",
#         legend.spacing.y = unit(0.0001, 'cm'),
#         legend.key.height = unit(0.01, "cm"), 
#         legend.key.width = unit(0.5, "cm"),
#         legend.box.just = "right",
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(fill = NA),
#         plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
#   # legend.spacing.x = unit(0.02, "cm"),
#   # legend.key.size = unit(0.3, 'cm')) +
#   guides(
#     fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
#     size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
#   labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count") 
# plot_GO_MF

plot_GO_MF <- ggplot(GO_list_MF) + 
  geom_vline(xintercept = -log10(0.05), color='red', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_MF$plotpoint, labels = GO_list_MF$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "skyblue", high = "#0719BC", breaks = c(round(min(GO_list_MF$Count, na.rm = TRUE)), round((max(GO_list_MF$Count, na.rm = TRUE) + min(GO_list_MF$Count, na.rm = TRUE) ) / 2), round(max(GO_list_MF$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(GO_list_MF$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=12, color='black'),
        axis.text.y = element_text(size=12, color='black'),
        axis.title = element_text(size=12, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        # plot.title = element_blank(),
        legend.title = element_text(size=10, color='black', hjust = 1),
        legend.text = element_text(size=8, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.195),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.01, 'cm'),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 1.5),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDR arms", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF

# GOMF_tableSave <- enGO_HDR_merged_MF %>%
#   dplyr::filter(p.adjust < 0.05) %>%
#   dplyr::arrange(desc(p.adjust)) %>% 
#   dplyr::select(-plotpoint) %>%
#   dplyr::mutate(geneID = gsub("/", ", ", geneID))


# Need to revert-back any simplifed Descritpions
# write.table(GOMF_tableSave, file = "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOMF_HDRgenes.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


 
# final_plot <- cowplot::plot_grid(
#   plot_GO_BP, plot_GO_MF,
#   ncol = 1,
#   align = "v",
#   axis = "lr",
#   labels = c("A", "B"),
#   label_size = 14,
#   label_fontface = "plain")
# final_plot


# ggsave("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOBP_GOMF_Cbrig.png", final_plot,  width = 7.5, height = 6, dpi = 600)



########### No enrichment in non-HDR arm domains ############
### BP ###
enGO_nHDR_merged <- clusterProfiler::enricher(
  gene = nHD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,NIC58),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = NIC58_merged_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_nHDR_merged)

dotplot(enGO_nHDR_merged, showCategory = 40, title = "BP HDRs")

enGO_nHDR_dt <- as.data.table(enGO_nHDR_merged@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())


insert_linebreaks <- function(s, width = 40) {
  sapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

GO_list_BP_nHDR <- enGO_nHDR_dt %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>% 
  # dplyr::slice_head(n=30) %>%
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number())

plot_GO_BP_nHDR <- ggplot(GO_list_BP_nHDR) +
  geom_vline(xintercept = -log10(0.05), color='red', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_BP_nHDR$plotpoint, labels = GO_list_BP_nHDR$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "skyblue", high = "#0719BC", breaks = c(round(min(GO_list_BP_nHDR$Count, na.rm = TRUE)), round((max(GO_list_BP_nHDR$Count, na.rm = TRUE) + min(GO_list_BP_nHDR$Count, na.rm = TRUE) ) / 2), round(max(GO_list_BP_nHDR$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(GO_list_BP_nHDR$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=12, color='black'),
        axis.text.y = element_text(size=12, color='black'),
        axis.title = element_text(size=12, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        # plot.title = element_blank(),
        legend.title = element_text(size=10, color='black', hjust = 1),
        legend.text = element_text(size=8, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.195),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.01, 'cm'),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 1.5),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in nHDR arms", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP_nHDR




### MF ###
enGO_nHDR_merged_MF <- clusterProfiler::enricher(
  gene = nHD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,NIC58),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = NIC58_merged_bgd_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_nHDR_merged_MF)

dotplot(enGO_nHDR_merged_MF, showCategory = 40, title = "MF nHDRs")

enGO_nHDR_merged_MF <- as.data.table(enGO_nHDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

insert_linebreaks <- function(s, width = 40) {
  sapply(s, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

GO_list_MF_nHDR <- enGO_nHDR_merged_MF %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>% 
  # dplyr::slice_head(n=30) %>%
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number()) 

plot_GO_MF_nHDR <- ggplot(GO_list_MF_nHDR) + 
  geom_vline(xintercept = -log10(0.05), color='red', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_MF_nHDR$plotpoint, labels = GO_list_MF_nHDR$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "skyblue", high = "#0719BC", breaks = c(round(min(GO_list_MF_nHDR$Count, na.rm = TRUE)), round((max(GO_list_MF_nHDR$Count, na.rm = TRUE) + min(GO_list_MF_nHDR$Count, na.rm = TRUE) ) / 2), round(max(GO_list_MF_nHDR$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(1, 10), name = "Fold enrichment", breaks = pretty(GO_list_MF_nHDR$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=12, color='black'),
        axis.text.y = element_text(size=12, color='black'),
        axis.title = element_text(size=12, color='black', face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        # plot.title = element_blank(),
        legend.title = element_text(size=10, color='black', hjust = 1),
        legend.text = element_text(size=8, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.195),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.01, 'cm'),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 1.5),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in nHDR arms", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF_nHDR
