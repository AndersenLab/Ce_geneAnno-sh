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


ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/N0_genes.tsv") %>%
  dplyr::rename(QX1410 = QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.protein, AF16 = c_briggsae.PRJNA10731.WS276.csq.ONLYPC.longest.protein, N2 = c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein)

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


#### Plotting ####
private_freq = (1/(length(strainCol_c2_u)))

classification <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq > private_freq & freq < 1 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>%
  dplyr::count(freq, class) %>%
  dplyr::mutate(percent = (n / sum(n)) * 100) 


gs_allOrtho <- ggplot(data = classification, aes(x = freq * 100, y = percent, fill = class)) + 
  geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
  scale_fill_manual(values = c(
    "core" = "green4",
    "accessory" = "#DB6333",
    "private" = "magenta3"
  ), 
  limits = c("core", "accessory", "private"),  # Manually ordering legend items
  guide = guide_legend(title = NULL) 
  ) +
  ylab("HOGs") + # longest isoform
  xlab("Frequency") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, color = 'black', face = 'bold'),
    legend.position = c(0.85, 0.8),
    plot.margin = margin(l = 20),
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=16, color = 'black'),
    axis.text = element_text(size=14, color = 'black')
  )
gs_allOrtho

HOG_class_count <- classification %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(n_HOG = sum(n)) %>%
  dplyr::ungroup()





# ======================================================================================================================================================================================== #
# GO - extract one/many-to-ones between QX and N2 and lift over GO terms #
# ======================================================================================================================================================================================== #

# Extracting QX - N2 one/many-to-ones
qn <- all_relations %>%
  dplyr::select(-AF16_count) %>%
  dplyr::filter(QX1410_count >= 1) %>%
  dplyr::filter(N2_count <= 1) 

oneoneqn <- ortho_genes_dd %>%
  dplyr::filter(HOG %in% qn$HOG) %>%
  dplyr::select(QX1410, N2)


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

N2_euk_bgd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output/euk_annotations/N2_longestIso_background_euk_GOannotations.tsv", col_names = c("transcript","GO", "term_name")) %>%
  dplyr::left_join(tran_gene, by = "transcript") %>%
  tidyr::separate_rows(GO, sep = ",") 

n2_gene_database <- gene2go %>%
  dplyr::select(WORMBASE, GO, ONTOLOGY) %>%
  dplyr::mutate(n2_gene = ifelse(WORMBASE %in% n2_wbGenes$gene, TRUE, FALSE)) %>%
  dplyr::mutate(n2_gene = ifelse(n2_gene == TRUE, WORMBASE, n2_gene)) %>%
  dplyr::mutate(eggNOG_GO = ifelse(n2_gene %in% N2_euk_bgd$gene, N2_euk_bgd$GO, NA))

merged  <- n2_gene_database %>%
  dplyr::select(n2_gene,GO,eggNOG_GO) %>%
  dplyr::filter(n2_gene != FALSE) %>% 
  dplyr::mutate(final_GO = dplyr::coalesce(GO, eggNOG_GO)) %>%
  dplyr::select(n2_gene,final_GO) %>%
  dplyr::filter(!is.na(final_GO)) %>%
  dplyr::rename(gene = n2_gene, GO = final_GO) # 14,345 N2 genes


# Lift over GO terms from N2 to QX1410
GO_brig <- oneoneqn %>%
  tidyr::separate_rows(QX1410, sep = ",\\s*") %>% 
  dplyr::left_join(merged, by = c("N2" = "gene"), relationship = "many-to-many") %>%
  dplyr::filter(!is.na(GO)) # 10,913


# ======================================================================================================================================================================================== #
# Lift over GO IDs from RBBH output #
# ======================================================================================================================================================================================== #
rbbh <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/RBBH/N2_QX1410_RBBH.tsv", col_names = c("N2","QX_rbbh"))
sumdf <- ortho_genes_dd %>%
  dplyr::select(HOG,QX1410, N2) %>%
  dplyr::left_join(GO_brig, by = "N2") %>%
  dplyr::select(-QX1410.y) %>%
  dplyr::rename(QX1410 = QX1410.x) %>%
  tidyr::separate_rows(N2, sep = ",") %>%
  dplyr::left_join(rbbh, by = "N2") 

rbbh_GO <- sumdf %>%
  dplyr::filter(is.na(GO) & !is.na(QX_rbbh)) %>%
  dplyr::select(HOG,QX1410,QX_rbbh,GO,N2) %>%
  dplyr::left_join(merged, by = c("N2" = "gene")) %>%
  dplyr::rename(GO_rbbh = GO.y)

join_final <- sumdf %>%
  dplyr::left_join(rbbh_GO, by = "N2") %>%
  dplyr::select(QX1410.x, N2, GO, QX_rbbh.x, GO_rbbh) %>%
  dplyr::rename(QX1410 = QX1410.x, QX_rbbh = QX_rbbh.x) %>%
  dplyr::filter(!is.na(QX1410) & !is.na(N2)) %>%
  tidyr::separate_rows(QX1410, sep = ",\\s") %>%
  dplyr::mutate(GO_final = ifelse(!is.na(GO), GO, GO_rbbh)) %>%
  dplyr::select(QX1410, N2, GO_final) %>%
  dplyr::filter(!is.na(GO_final)) # 11,516 QX1410 genes


# ======================================================================================================================================================================================== #
# After adding output from RBBH, then look at QX1410 - AF16 one/many-to-ones and then convert AF16 gene name to N2 WBGeneID and pull any GO IDs #
# ======================================================================================================================================================================================== #
cb_alias <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/CB_N2_alias.tsv", col_names = c("AF_WBGeneID","alias"))
N2_alias <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/N2_aliasNames.tsv", col_names = c("N2_WBGeneID","alias")) %>%
  dplyr::select(alias,N2_WBGeneID)


# Isolating QX1410 - AF16 one/many-to-one
cbcb <- all_relations %>%
  dplyr::select(-N2_count) %>%
  dplyr::filter(QX1410_count >= 1) %>%
  dplyr::filter(AF16_count <= 1) 

# Adding AF16 alias names and left-joining with N2 alias' for GO ID lift over
cbcb_genes <- ortho_genes_dd %>%
  dplyr::filter(HOG %in% cbcb$HOG) %>%
  dplyr::select(QX1410, AF16) %>%
  dplyr::left_join(cb_alias, by = c("AF16" = "AF_WBGeneID")) %>%
  dplyr::left_join(N2_alias, by = "alias") %>%
  tidyr::separate_rows(QX1410, sep = ',\\s') %>%
  dplyr::left_join(join_final, by = "QX1410") %>%
  dplyr::filter(is.na(GO_final) & is.na(N2) & !is.na(N2_WBGeneID)) %>%
  dplyr::left_join(merged, by = c("N2_WBGeneID" = "gene")) %>%
  dplyr::rename(GO_fromAlias = GO)

final_GO <- ortho_genes_dd %>%
  dplyr::select(QX1410) %>%
  tidyr::separate_rows(QX1410, sep = ',\\s') %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::left_join(join_final, by = "QX1410") %>%
  dplyr::left_join(cbcb_genes, by = "QX1410") %>%
  dplyr::select(QX1410, N2.x, GO_final.x, N2_WBGeneID, GO_fromAlias) %>%
  dplyr::mutate(N2 = ifelse(!is.na(N2.x), N2.x, N2_WBGeneID)) %>%
  dplyr::mutate(GO = ifelse(!is.na(GO_final.x), GO_final.x, GO_fromAlias)) %>%
  dplyr::select(QX1410, N2, GO) %>%
  dplyr::filter(!is.na(N2))





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


nHDR_armDomain <- readr::read_csv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/briggsae/multi_reference/2024_release/QX1410/chromosome_domain_Cbriggsae.csv") %>%
  dplyr::rename(Chromosome=chrom) %>%
  dplyr::select(Chromosome,left,right) %>%
  dplyr::mutate(left=left*1e3,right=right*1e3)

hdr_regions <- readr::read_tsv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/briggsae/multi_reference/reference_realign/HDR_CB_allStrain_5kbclust_20250613.tsv") 

# hdr_filt <- hdr_regions %>% 
  # dplyr::filter(source == "QX1410")

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

# Read in GFF for extracting QX1410 genes in nHDRs and HDRs
gffCB <-  ape::read.gff("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.gff3")

gffCB_genes <- gffCB %>%
  dplyr::filter(type=="gene") %>%
  dplyr::filter(grepl("biotype=protein_coding", attributes)) %>%
  tidyr::separate(attributes, into=c("ID","rest"), sep = "ID=gene:") %>%
  tidyr::separate(rest, into=c("keep","nope"), sep = ";biotype=") %>%
  dplyr::select(-ID,-source,-type,-score,-strand,-phase,-nope) %>%
  dplyr::rename(CHROM = seqid, QX1410 = keep) %>%
  dplyr::left_join(final_GO, by = "QX1410")

gffFinal <- gffCB_genes %>%
  dplyr::filter(!is.na(N2)) %>%
  dplyr::select(-GO) %>%
  dplyr::distinct(QX1410, .keep_all = T)


######### THIS IS DONE WITH N2 GENES AND N2 BACKGROUND ############
# HDreg_arm <- ldply(reglist,data.frame)
# HDreg_center <- ldply(reglist_center,data.frame)
# HDgene <- getGenesFromBed(gffFinal,HDreg_arm)
# HD_gene_vector <- HDgene[[]]
# HD_n2_genes <- (do.call(rbind, HDgene))
# # testHDR_genes <- dplyr::select(HD_n2_genes, N2) 
# HD_gene_vector <- unique(HD_n2_genes$N2) # pulling a vector list of N2 genes that are in HDRs on QX1410 chromosomal arms
# 
# nHDreg <- classifyingArms(HDreg_arm %>% dplyr::rename(minStart=start,maxEnd=end), nHDR_armDomain) %>% dplyr::select(-domain)
# nHDgene <- getGenesFromBed(gffFinal,nHDreg)
# nHD_n2_genes <- do.call(rbind, nHDgene)
# # testnHDR_genes <- dplyr::select(nHD_n2_genes, N2) 
# nHD_gene_vector <- unique(nHD_n2_genes$N2) # pulling a vector list of N2 genes that are NOT in HDRs on QX1410 chromosomal arms
# 
# 
# posplot <- ggplot() + 
#   geom_rect(data=centers %>% dplyr::rename(CHROM=Chromosome), aes(xmin=center_start,xmax=center_end,ymin=-1,ymax=1,fill="nHDR_center")) +
#   geom_rect(data=rbind(HDreg_arm %>% dplyr::mutate(class="HDR_arm"),nHDreg %>% dplyr::mutate(class="nHDR_arm"), HDreg_center %>% dplyr::mutate(class="HDR_center")), aes(xmin=start,xmax=end,ymin=-1,ymax=1,fill=class)) +
#   facet_wrap(~CHROM,ncol=1,scales="free_x") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid = element_blank()) +
#   scale_fill_manual(values=c("nHDR_center"="lightblue","nHDR_arm"="deepskyblue3","HDR_center"="pink","HDR_arm"="firebrick3"))+
#   xlab("Physical postion (Mb)")
# posplot
# 
# posplot2 <- ggplot() + 
#   geom_rect(data=HD_n2_genes, aes(xmin=start,xmax=end,ymin=-1,ymax=1)) +
#   facet_wrap(~CHROM,ncol=1,scales="free_x") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid = element_blank()) +
#   xlab("Physical postion (Mb)")
# posplot2
# 
# 
# 
# N2_merged_bgd_genes <- unique(merged$gene) # a vector of genes
# GO_terms_merged <- AnnotationDbi::select(GO.db, 
#                                          keys=unique(merged$GO),   
#                                          columns = c("TERM", "DEFINITION", "ONTOLOGY"),
#                                          keytype="GOID") %>% 
#   dplyr::rename(TERM = GOID, TERM_NAME = TERM)
# 
# merged_ont <- merged %>%
#   dplyr::left_join(GO_terms_merged, by = c("GO" = "TERM")) %>%
#   dplyr::mutate(inHDR = ifelse(gene %in% testHDR_genes$N2, TRUE, FALSE))
# 
# 
# # HDR genes in arms against all N2 background set
# enGO_HDR_merged <- clusterProfiler::enricher(
#   gene = HD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_merged_bgd_genes,
#   pvalueCutoff = 0.5,
#   # pAdjustMethod = "BH",
#   qvalueCutoff = 0.5,
# )
# 
# head(enGO_HDR_merged)
# 
# dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")
# 
# # enGO_HDR_genes <- clusterProfiler::setReadable(enGO_HDR_merged, keyType = "WORMBASE", OrgDb = org.Ce.eg.db)
# 
# 
# # nHDRs genes in arms against all N2 background set
# enGO_nHDR_merged <- clusterProfiler::enricher(
#   gene = nHD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "CC") %>% dplyr::select(GO,gene),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "CC") %>% dplyr::select(TERM,TERM_NAME),
#   universe = N2_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "bonferroni",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_nHDR_merged)
# 
# dotplot(enGO_nHDR_merged, showCategory = 40, title = "BP QX1410 nHDRs")



















######### THIS IS DONE WITH QX1410 GENES AND QX1410 BACKGROUND ############
HDreg_arm <- ldply(reglist,data.frame)
HDreg_center <- ldply(reglist_center,data.frame)
HDgene <- getGenesFromBed(gffFinal,HDreg_arm)
HD_gene_vector <- HDgene[[]]
HD_QX_genes <- (do.call(rbind, HDgene))
# testHDR_genes <- dplyr::select(HD_n2_genes, N2) 
HD_gene_vector <- unique(HD_QX_genes$QX1410) # pulling a vector list of N2 genes that are in HDRs on QX1410 chromosomal arms

nHDreg <- classifyingArms(HDreg_arm %>% dplyr::rename(minStart=start,maxEnd=end), nHDR_armDomain) %>% dplyr::select(-domain)
nHDgene <- getGenesFromBed(gffFinal,nHDreg)
nHD_QX_genes <- do.call(rbind, nHDgene)
# testnHDR_genes <- dplyr::select(nHD_n2_genes, N2) 
nHD_gene_vector <- unique(nHD_QX_genes$QX1410) # pulling a vector list of N2 genes that are NOT in HDRs on QX1410 chromosomal arms


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
  geom_rect(data=HD_QX_genes, aes(xmin=start,xmax=end,ymin=-1,ymax=1)) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  xlab("Physical postion (Mb)")
posplot2



QX1410_merged_bgd_genes <- unique(final_GO$QX1410) #  12,100 genes
GO_terms_merged <- AnnotationDbi::select(GO.db, 
                                         keys=unique(final_GO$GO),   
                                         columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                         keytype="GOID") %>% 
  dplyr::rename(TERM = GOID, TERM_NAME = TERM) # 6,960

merged_ont <- final_GO %>%
  dplyr::left_join(GO_terms_merged, by = c("GO" = "TERM")) 

# HDR genes in arms against all N2 background set
enGO_HDR_merged <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = QX1410_merged_bgd_genes,
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
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number())

# GO_list_BP <- enGO_HDR_dt %>%
#   dplyr::filter(p.adjust < 0.05) %>%
#   dplyr::arrange(desc(p.adjust)) %>%  
#   dplyr::mutate(Description = factor(Description, levels = unique(Description)), plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(GO_list_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_BP$plotpoint, labels = GO_list_BP$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(GO_list_BP$Count, na.rm = TRUE)), round((max(GO_list_BP$Count, na.rm = TRUE) + min(GO_list_BP$Count, na.rm = TRUE) ) / 2), round(max(GO_list_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.1, 3), name = "Fold enrichment", breaks = pretty(GO_list_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=6, color='black'), 
        axis.text.y = element_text(size=6, color='black'),
        axis.title = element_text(size=7, color='black', face = 'bold'), 
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.title = element_text(size=4, color='black', hjust = 1), 
        legend.text = element_text(size=3, color='black'), 
        plot.title = element_blank(),
        legend.position = "inside", 
        legend.position.inside = c(0.8, 0.18),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),   
        legend.key.width = unit(0.5, "cm"),
        legend.justification = c("right", "top"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count") 
plot_GO_BP


GOBP_tableSave <- GO_list_BP %>%
  dplyr::select(-plotpoint) %>%
  dplyr::mutate(geneID = gsub("/", ", ", geneID))

# write.table(GOBP_tableSave, file = "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOBP_HDRgenes.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = QX1410_merged_bgd_genes,
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
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number()) 

GO_list_MF <- GO_list_MF %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen","oxidoreductase activity (1)", Description), 
                Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (2)", Description),
                Description = gsub("transmitter-gated ion channel activity involved in regulation of postsynaptic membrane potential", "transmitter-gated ion channel (1)", Description))


# oxidoreductase activity (1) = oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen
# oxidoreductase activity (2) = oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen
# transmitter-gated ion channel (1) = transmitter-gated ion channel activity involved in regulation of postsynaptic membrane potential



plot_GO_MF <- ggplot(GO_list_MF) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_MF$plotpoint, labels = GO_list_MF$Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(GO_list_MF$Count, na.rm = TRUE)), round((max(GO_list_MF$Count, na.rm = TRUE) + min(GO_list_MF$Count, na.rm = TRUE) ) / 2), round(max(GO_list_MF$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.1, 3), name = "Fold enrichment", breaks = pretty(GO_list_MF$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=6, color='black'), 
        axis.text.y = element_text(size=6, color='black'),
        axis.title = element_text(size=7, color='black', face = 'bold'), 
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.title = element_text(size=4, color='black', hjust = 1), 
        legend.text = element_text(size=3, color='black'), 
        plot.title = element_blank(),
        legend.position = "inside", 
        legend.position.inside = c(0.8, 0.18),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"), 
        legend.key.width = unit(0.5, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 20, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count") 
plot_GO_MF

GOMF_tableSave <- GO_list_MF %>%
  dplyr::select(-plotpoint, -Wrapped_Description) %>%
  dplyr::mutate(geneID = gsub("/", ", ", geneID))

# write.table(GOMF_tableSave, file = "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOMF_HDRgenes.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



final_plot <- cowplot::plot_grid(
  plot_GO_BP, plot_GO_MF, 
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("A", "B"),
  label_size = 14,       
  label_fontface = "plain")
final_plot


ggsave("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOBP_GOMF_Cbrig.png", final_plot,  width = 7.5, height = 6, dpi = 600)


