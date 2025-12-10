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
#library(broom) # added by Nic


ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/Orthogroups_N0_genes.tsv") %>%
  dplyr::rename(QX1410 = QX1410.longest.prot, AF16 = c_briggsae.PRJNA10731.WS276.csq.ONLYPC.longest.protein, N2 = c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.protein)

strainCol <- colnames(ortho_genes_dd)
colnames(ortho_genes_dd) <- strainCol

ortho_count <- ortho_genes_dd

strainCol_c2_u <- strainCol[!strainCol %in% c("Orthogroup")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations <- ortho_count %>%
  dplyr::select(Orthogroup, dplyr::contains("_count"))


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
  ylab("OGs") + # longest isoform
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

# Prepping for enrichment analysis

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

### Collapsing HDRs found among all strains ###
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
    # print(i)
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

ggplot(data = nHDR_armDomain) + 
  geom_rect(aes(xmin = left / 1e6, xmax = right / 1e6, ymin = 0.5, ymax = 1.5, fill = 'red')) +
  facet_wrap(~Chromosome, scales = 'free') + 
  theme(
    legend.position = 'none')

hdr_regions <- readr::read_tsv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/briggsae/multi_reference/reference_realign/HDR_CB_allStrain_5kbclust_20250930.tsv") %>%
  dplyr::filter(STRAIN != "ECA1605" & STRAIN != "ECA1559")

hdr_regions <- hdr_regions %>%
  dplyr::filter(source == "QX1410")
# "QX1410"   "BRC20530" "NIC1660"  "BRC20492" "QG2902"   "NIC1667"  "ED3102"   "ECA2670"  "ECA2666"  "JU1348"   "QG1005"   "QG2964"

HDreg_all <- ldply(getRegFreq(hdr_regions %>% dplyr::select(CHROM,minStart,maxEnd,STRAIN) %>% dplyr::arrange(CHROM,minStart) %>%dplyr::group_split(CHROM)), data.frame) %>% dplyr::select(-STRAIN) %>% dplyr::rename(start=minStart,end=maxEnd)

centers <- nHDR_armDomain %>%
  dplyr::group_by(Chromosome) %>%
  dplyr::mutate(center_end=lead(left)) %>%
  dplyr::filter(!is.na(center_end)) %>%
  dplyr::select(-left) %>%
  dplyr::rename(center_start=right)

# Filtering for HDRs NOT in chromosomal centers
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
HDreg_arm_and_tip <- ldply(reglist,data.frame)

# Filtering for HDRs in arm domains
reglist2 <- list()
for (i in 1:nrow(nHDR_armDomain)) {
  astart <- nHDR_armDomain[i,]$left
  aend <- nHDR_armDomain[i,]$right
  aend_lead <- n
  achrom <- nHDR_armDomain[i,]$Chromosome
  temp <- HDreg_arm_and_tip %>%
    dplyr::filter(CHROM==achrom) %>%
    dplyr::filter((start >= astart & end <= aend) | (start <= astart & end >= astart) | (start <= aend & end >= aend))
  
  reglist2[[i]] <- temp
}

# Read in GFF for extracting QX1410 genes in nHDRs and HDRs
gffCB <-  ape::read.gff("/vast/eande106/projects/Nicolas/CBCN_comp/tree/gene_level_20251006/updated_QX/c_briggsae.QX1410_20250929.csq.longest.gff")

gffFinal <- gffCB %>%
  dplyr::filter(type=="gene") %>%
  dplyr::filter(grepl("biotype=protein_coding", attributes)) %>%
  tidyr::separate(attributes, into=c("ID","rest"), sep = "ID=gene:") %>%
  tidyr::separate(rest, into=c("keep","nope"), sep = ";biotype=") %>%
  dplyr::select(-ID,-source,-type,-score,-strand,-phase,-nope) %>%
  dplyr::rename(CHROM = seqid, QX1410 = keep) %>%
  # dplyr::left_join(final_GO, by = "QX1410")
  dplyr::distinct(QX1410, .keep_all = T)

# gffFinal <- gffCB_genes %>%
#   dplyr::filter(!is.na(N2)) %>%
#   dplyr::select(-GO) %>%
#   dplyr::distinct(QX1410, .keep_all = T)

# HDR genes in arm domains
HDreg_arm <- ldply(reglist2, data.frame)
HDreg_center <- ldply(reglist_center,data.frame)
HDgene <- getGenesFromBed(gffFinal,HDreg_arm)
HD_gene_vector <- HDgene[[]]
HD_QX_genes <- (do.call(rbind, HDgene))
HD_gene_vector <- unique(HD_QX_genes$QX1410) 


## Extracting ALL HDR genes ##
# HD_all <- dplyr::bind_rows(HDreg_arm, HDreg_center)
# HDgene_all <- getGenesFromBed(gffFinal,HD_all)
# HD_all_gene_vector <- HDgene_all[[]]
# HD_QX_genes_all <- (do.call(rbind, HDgene_all))
# HD_all_gene_vector <- unique(HD_QX_genes_all$QX1410)

nHDreg <- classifyingArms(HDreg_arm_and_tip %>% dplyr::rename(minStart=start,maxEnd=end), nHDR_armDomain) %>% dplyr::select(-domain)
nHDgene <- getGenesFromBed(gffFinal,nHDreg)
nHD_QX_genes <- do.call(rbind, nHDgene)
nHD_gene_vector <- unique(nHD_QX_genes$QX1410) 





# ========================================================================================================================================================================================== #
# All HDRs - not JUST arms 
# hdr_regions_cleaned <- hdr_regions %>% dplyr::rename(start = minStart, end = maxEnd) %>% dplyr::select(CHROM, start, end)
# hdr_genes <- getGenesFromBed(gffFinal, hdr_regions_cleaned) # 
# all_HDR_genes <- do.call(rbind, hdr_genes)
# all_HDR_genes_vector <- unique(all_HDR_genes$QX1410)
#   
# all_non_HDR_genes
# 
# # ALL HDR GENES
# posplot <- ggplot() + 
#   geom_rect(data = hdr_regions, aes(xmin=minStart,xmax=maxEnd,ymin=-1,ymax=1), fill = 'black') +
#   # geom_rect(data = all_HDR_genes, aes(xmin = , xmax = , ymin = -1, ymax = 1), fill = 'red') +
#   facet_wrap(~CHROM,ncol=1,scales="free_x") +
#   theme_bw() +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid = element_blank()) +
#   scale_fill_manual(values=c("nHDR_center"="lightblue","nHDR_arm"="deepskyblue3","HDR_center"="pink","HDR_arm"="firebrick3"))+
#   xlab("Physical postion (Mb)")
# posplot
# ========================================================================================================================================================================================== #




posplot <- ggplot() + 
  geom_rect(data=centers %>% dplyr::rename(CHROM=Chromosome), aes(xmin=center_start,xmax=center_end,ymin=-1,ymax=1,fill="nHDR_center")) +
  # geom_rect(data = hdr_regions, aes(xmin=minStart,xmax=maxEnd,ymin=-1,ymax=1)) +
  geom_rect(data=rbind(HDreg_arm %>% dplyr::mutate(class="HDR_arm"),nHDreg %>% dplyr::mutate(class="nHDR_arm"), HDreg_center %>% dplyr::mutate(class="HDR_center")), aes(xmin=start,xmax=end,ymin=-1,ymax=1,fill=class)) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values=c("nHDR_center"="lightblue","nHDR_arm"="deepskyblue3","HDR_center"="pink","HDR_arm"="firebrick3"))+
  xlab("Physical postion (Mb)")
posplot

posplot1 <- ggplot() + 
  geom_rect(data=centers %>% dplyr::rename(CHROM=Chromosome), aes(xmin=center_start,xmax=center_end,ymin=-2,ymax=2, fill="center")) +
  geom_rect(data = nHDR_armDomain %>% dplyr::rename(CHROM=Chromosome), aes(xmin=left,xmax=right,ymin=-2,ymax=2, fill = "arm_region")) +
  geom_rect(data = HDreg_arm, aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = "HDR_arm_region")) +
  geom_rect(data = nHDreg, aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = "nHDR_arm_region")) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  # scale_fill_manual(values=c("nHDR_center"="lightblue","nHDR_arm"="deepskyblue3","HDR_center"="pink","HDR_arm"="firebrick3"))+
  xlab("Physical postion (Mb)")
posplot1

posplot2 <- ggplot() + 
  geom_rect(data=HD_QX_genes, aes(xmin=start,xmax=end,ymin=-1,ymax=1)) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  xlab("Physical postion (Mb)")
posplot2


# Genes in arm domains - serve as background set
arms <- nHDR_armDomain %>% dplyr::rename(CHROM=Chromosome,start=left,end=right)
all_arm_genes <- getGenesFromBed(gffFinal,arms)
all_arm_genes_df <- as.data.frame((do.call(rbind, all_arm_genes)))
arm_genes <- all_arm_genes_df$QX1410

ggplot(data = nHDR_armDomain) + 
  geom_rect(aes(xmin = left / 1e6, xmax = right / 1e6, ymin = 0.5, ymax = 1.5), fill = 'black') +
  geom_rect(data = all_arm_genes_df %>% dplyr::rename(Chromosome = CHROM), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.5, ymax = 1.5), fill = 'gold') +
  # geom_rect(data = HD_QX_genes %>% dplyr::rename(Chromosome = CHROM), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.75, ymax = 1.25), fill = "blue") +
  facet_wrap(~Chromosome, scales = 'free') + 
  theme(
    legend.position = 'none')

ggplot(data = nHDR_armDomain) + 
  geom_rect(aes(xmin = left / 1e6, xmax = right / 1e6, ymin = 0.5, ymax = 1.5), fill = 'black') +
  geom_rect(data = all_arm_genes_df %>% dplyr::rename(Chromosome = CHROM), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.5, ymax = 1.5), fill = 'gold') +
  # geom_rect(data = HD_QX_genes %>% dplyr::rename(Chromosome = CHROM), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.75, ymax = 1.25), fill = "blue") +
  facet_wrap(~Chromosome, scales = 'free') + 
  theme(
    legend.position = 'none')




#==============================================================================================================================================================================================================================#

# INTERPROSCAN - all genes as background and ALL HDR genes are enrichment set

#==============================================================================================================================================================================================================================#
# Run InterProScan:
# interproscan.sh \
# --formats TSV \
# --input /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/prep/QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.protein.fa \
# --goterms \
# --aplications Pfam,SMART,TIGRFAM,SUPERFAMILY,CDD,Gene3D,FunFam,PANTHER,PIRSF,PIRSR,ProSiteProfiles,ProSitePatterns,SFLD,Hamap,Coils,AntiFam \
# --cpu 48 \
# --iprlookup \
# --disable-precalc \
# --output-file-base $output/QX1410_InterProScan \
# --tempdir $TMP

# awk -F'\t' -v OFS='\t' '{print $1,$5,$6,$12,$13,$14}' QX1410_InterProScan.tsv > QX1410_IPR_GO.tsv

allQX <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/QX1410_tran_gene_updatedGFF_20251007.tsv", col_names = c("tran","QX1410")) %>%
  dplyr::mutate(tran = paste0("transcript:",tran))


### Analysis of all HDR genes against ALL background genes, and analysis of all nHDR genes against ALL background genes
ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/output/QX_IPR_allApps_20251006.tsv", col_names = c("tran", "MD5_digest", "seq_length", "app", "signature_accession", "signature_description", "start", "end", "score", "status", "date", "IPR_accession","IPR_description","GO", "pathways")) %>%
  dplyr::left_join(allQX, by = 'tran') %>%
  dplyr::select(-tran) %>%
  dplyr::select(QX1410,signature_accession,signature_description,IPR_accession,IPR_description,GO) 

# ipr_table <- ipr %>%
#   dplyr::select(QX1410, IPR_accession, IPR_description, GO) %>%
#   dplyr::filter(IPR_accession != '-' | GO != "-") %>%
#   dplyr::rename(GO_ID = GO)
# 
# write.table(ipr_table, "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/output/QX_IPR_cleaned_updated_20251007.tsv", sep = '\t', quote = F, col.names = T, row.names = F)



# # Unique InterPro per gene
# ipr_gene <- ipr %>%
#   dplyr::filter(!is.na(IPR_description) & IPR_description != "-") %>%
#   dplyr::select(QX1410, IPR_accession, IPR_description) %>%
#   dplyr::distinct(QX1410,IPR_accession, IPR_description) # 15,289 genes!
# 
# # Define universe & HDR membership (annotated-only universe) 
# univ_genes <- unique(ipr_gene$QX1410)
# hdr_genes  <- intersect(all_HDR_genes_vector, univ_genes) # 3,427 genes 
# 
# N <- length(univ_genes)
# n <- length(hdr_genes)
# 
# # Counts per IPR (k = in universe, x = in HDR subset) 
# k_tbl <- ipr_gene %>%
#   dplyr::count(IPR_accession, name = "k")
# 
# x_tbl <- ipr_gene %>%
#   dplyr::filter(QX1410 %in% hdr_genes) %>%
#   dplyr::count(IPR_accession, name = "x")
# 
# desc_tbl <- ipr_gene %>%
#   dplyr::distinct(IPR_accession, IPR_description)
# 
# # Hypergeometric enrichment (one-sided)
# ipr_enrichment <- k_tbl %>%
#   dplyr::left_join(x_tbl, by = "IPR_accession") %>%
#   dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
#   dplyr::mutate(
#     pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
#     expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
#     enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
#     # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
#     OR = {
#       a <- x + 0.5                                # HDR & has IPR
#       b <- (n - x) + 0.5                          # HDR & no IPR
#       c <- (k - x) + 0.5                          # non-HDR & has IPR
#       d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
#       (a / b) / (c / d)
#     },
#     FDR_p.adjust = stats::p.adjust(pval, method = "BH")
#   ) %>%
#   dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
#   dplyr::mutate(N = N, n = n) %>%
#   dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
#   dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))
# 
# ipr_sig <- ipr_enrichment %>%
#   dplyr::filter(FDR_p.adjust < 0.05)
# 
# ipr_sig %>% dplyr::slice_head(n = 20)
# 
# ipr_sig_gene_collapsed <- ipr_gene %>%
#   dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
#   dplyr::group_by(IPR_accession) %>%
#   dplyr::summarise(
#     IPR_description = dplyr::first(stats::na.omit(IPR_description)),
#     n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
#     genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
#     n_genes_all = dplyr::n_distinct(QX1410),
#     genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
#     .groups = "drop") %>%
#   dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
#   dplyr::mutate(Region = "hyper-divergent regions")
# 
# 
# # Now for non-HDR arm genes
# univ_genes2 <- unique(ipr_gene$QX1410)
# nhdr_genes  <- intersect(all_non_HDR_genes, univ_genes2) # 3,427 genes 
# 
# N <- length(univ_genes2)
# n <- length(nhdr_genes)
# 
# # Counts per IPR (k = in universe, x = in HDR subset) 
# k_tbl2 <- ipr_gene %>%
#   dplyr::count(IPR_accession, name = "k")
# 
# x_tbl2 <- ipr_gene %>%
#   dplyr::filter(QX1410 %in% nhdr_genes) %>%
#   dplyr::count(IPR_accession, name = "x")
# 
# desc_tbl2 <- ipr_gene %>%
#   dplyr::distinct(IPR_accession, IPR_description)
# 
# # Hypergeometric enrichment (one-sided)
# ipr_enrichment_nHDR <- k_tbl2 %>%
#   dplyr::left_join(x_tbl2, by = "IPR_accession") %>%
#   dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
#   dplyr::mutate(
#     pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
#     expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
#     enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
#     # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
#     OR = {
#       a <- x + 0.5                                # HDR & has IPR
#       b <- (n - x) + 0.5                          # HDR & no IPR
#       c <- (k - x) + 0.5                          # non-HDR & has IPR
#       d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
#       (a / b) / (c / d)
#     },
#     FDR_p.adjust = stats::p.adjust(pval, method = "BH")
#   ) %>%
#   dplyr::left_join(desc_tbl2, by = "IPR_accession") %>%
#   dplyr::mutate(N = N, n = n) %>%
#   dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
#   dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))
# 
# ipr_sig2 <- ipr_enrichment_nHDR %>%
#   dplyr::filter(FDR_p.adjust < 0.05)
# 
# ipr_sig2 %>% dplyr::slice_head(n = 20)
# 
# ipr_sig_gene_collapsed2 <- ipr_gene %>%
#   dplyr::filter(IPR_accession %in% ipr_sig2$IPR_accession) %>%
#   dplyr::group_by(IPR_accession) %>%
#   dplyr::summarise(
#     IPR_description = dplyr::first(stats::na.omit(IPR_description)),
#     n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
#     genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
#     n_genes_all = dplyr::n_distinct(QX1410),
#     genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
#     .groups = "drop") %>%
#   dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
#   dplyr::mutate(Region = "non HDRs")
# 
# 
# binded <- ipr_sig_gene_collapsed %>% dplyr::bind_rows(ipr_sig_gene_collapsed2) %>% dplyr::arrange(FDR_p.adjust) 
# 
# data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())
# 
# plot_ipr <- ggplot(data_plt) +
#   geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
#   geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
#   scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
#   scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
#   scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
#   scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
#   theme(axis.text.x = element_text(size=8, color='black'),
#         axis.text.y = element_text(size=8, color='black'),
#         axis.title = element_text(size=9, color='black', face = 'bold'),
#         axis.title.x = element_blank(),
#         # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#         plot.title = element_blank(),
#         legend.title = element_text(size=6.5, color='black', hjust = 1),
#         legend.text = element_text(size=5.5, color='black', hjust = 1),
#         legend.position = "inside",
#         legend.position.inside = c(0.76, 0.19),
#         legend.direction = "horizontal", legend.box = "vertical",
#         legend.spacing.y = unit(0.0001, 'cm'),
#         legend.key.height = unit(0.01, "cm"),
#         legend.key.width = unit(0.5, "cm"),
#         legend.box.just = "right",
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(fill = NA),
#         plot.margin = margin(b = 5, t = 10, r = 10, l = 22, unit = "pt")) +
#   # legend.spacing.x = unit(0.02, "cm"),
#   # legend.key.size = unit(0.3, 'cm')) +
#   guides(
#     fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
#     size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
#     shape = 'none') +
#   # shape = guide_legend(nrow=1, order = 3, title.position = "top", title.hjust = 1, force = TRUE)) +
#   labs(title = "Enriched IPR terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
# plot_ipr







#==============================================================================================================================================================================================================================#

# INTERPROSCAN - just ARM genes as background

#==============================================================================================================================================================================================================================#
ipr_gene <- ipr %>%
  dplyr::filter(!is.na(IPR_description) & IPR_description != "-") %>%
  dplyr::select(QX1410, IPR_accession, IPR_description) %>%
  dplyr::distinct(QX1410,IPR_accession, IPR_description) %>% # 15,289
  dplyr::filter(QX1410 %in% arm_genes) # 5,301 genes!

# Define universe & HDR membership (annotated-only universe) 
univ_genes <- unique(ipr_gene$QX1410)
hdr_genes  <- intersect(HD_gene_vector, univ_genes) # 3,087 genes 

N <- length(univ_genes)
n <- length(hdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl <- ipr_gene %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl <- ipr_gene %>%
  dplyr::filter(QX1410 %in% hdr_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl <- ipr_gene %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment <- k_tbl %>%
  dplyr::left_join(x_tbl, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig <- ipr_enrichment %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed <- ipr_gene %>%
  dplyr::filter(IPR_accession %in% ipr_sig$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
    genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(QX1410),
    genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "hyper-divergent regions")


# Now for non-HDR arm genes
univ_genes2 <- unique(ipr_gene$QX1410)
nhdr_genes  <- intersect(nHD_gene_vector, univ_genes2) # 3,427 genes 

N <- length(univ_genes2)
n <- length(nhdr_genes)

# Counts per IPR (k = in universe, x = in HDR subset) 
k_tbl2 <- ipr_gene %>%
  dplyr::count(IPR_accession, name = "k")

x_tbl2 <- ipr_gene %>%
  dplyr::filter(QX1410 %in% nhdr_genes) %>%
  dplyr::count(IPR_accession, name = "x")

desc_tbl2 <- ipr_gene %>%
  dplyr::distinct(IPR_accession, IPR_description)

# Hypergeometric enrichment (one-sided)
ipr_enrichment_nHDR <- k_tbl2 %>%
  dplyr::left_join(x_tbl2, by = "IPR_accession") %>%
  dplyr::mutate(x = tidyr::replace_na(x, 0L)) %>%
  dplyr::mutate(
    pval = stats::phyper(q = x - 1, m = k, n = N - k, k = n, lower.tail = FALSE),
    expected = (n * k) / N, # if IPR genes are randomly distributed, you’d expect this many HDR genes to carry the IPR.
    enrich_ratio = dplyr::if_else(expected > 0, x / expected, NA_real_), # (x HDR with IPR / k background with IPR) / (n HDR genes / N background genes)
    # odds ratio with Haldane–Anscombe correction (adding 0.5 to each cell to avoid infinities)
    OR = {
      a <- x + 0.5                                # HDR & has IPR
      b <- (n - x) + 0.5                          # HDR & no IPR
      c <- (k - x) + 0.5                          # non-HDR & has IPR
      d <- (N - n - (k - x)) + 0.5                # non-HDR & no IPR
      (a / b) / (c / d)
    },
    FDR_p.adjust = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(desc_tbl2, by = "IPR_accession") %>%
  dplyr::mutate(N = N, n = n) %>%
  dplyr::select(IPR_accession, IPR_description, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust) %>%
  dplyr::arrange(FDR_p.adjust, dplyr::desc(enrich_ratio))

ipr_sig2 <- ipr_enrichment_nHDR %>%
  dplyr::filter(FDR_p.adjust < 0.05)

ipr_sig2 %>% dplyr::slice_head(n = 20)

ipr_sig_gene_collapsed2 <- ipr_gene %>%
  dplyr::filter(IPR_accession %in% ipr_sig2$IPR_accession) %>%
  dplyr::group_by(IPR_accession) %>%
  dplyr::summarise(
    IPR_description = dplyr::first(stats::na.omit(IPR_description)),
    n_genes_HDR = dplyr::n_distinct(QX1410[QX1410 %in% hdr_genes]),
    genes_HDR   = paste(sort(unique(QX1410[QX1410 %in% hdr_genes])), collapse = ", "),
    n_genes_all = dplyr::n_distinct(QX1410),
    genes_all   = paste(sort(unique(QX1410)), collapse = ", "),
    .groups = "drop") %>%
  dplyr::left_join(ipr_sig2 %>% dplyr::select(IPR_accession, x, k, n, N, expected, enrich_ratio, OR, pval, FDR_p.adjust), by = "IPR_accession") %>%
  dplyr::mutate(Region = "non HDRs")


binded <- ipr_sig_gene_collapsed %>% dplyr::bind_rows(ipr_sig_gene_collapsed2) %>% dplyr::arrange(FDR_p.adjust) 

data_plt <- binded %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) + #### CHANGE TO SHAPE = REGION in aes FOR HDR AND nHDR IN SAME PLOT
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_shape_manual(values = c("hyper-divergent regions" = 21, "non HDRs" = 22)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size = 6, color='black', hjust = 1),
        legend.text = element_text(size = 5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.2),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, t = 10, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill = guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE),
    shape = 'none') +
  # shape = guide_legend(nrow=1, order = 3, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/final_plotsAndData/IPR_armBackground_HDRandnHDR_Cbrig.png", plot_ipr,  width = 7.5, height = 4.5, dpi = 600)



#==============================================================================================================================================================================================================================#

# INTERPROSCAN (Gene Ontology) - just ARM genes as background - HDR arm domain genes

#==============================================================================================================================================================================================================================#
go_ipr <- ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(QX1410, GO) %>% # 11,756 genes
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())


### Now with only arms as the background, not the entire genome
go_ipr_arms <- go_ipr %>% dplyr::filter(QX1410 %in% arm_genes)
IPR_GO_bckgrd_arms <- unique(go_ipr_arms$QX1410) # 3,905 genes


how_many_HDR_GO_arm_genes <- go_ipr_arms %>% dplyr::filter(QX1410 %in% HD_gene_vector) # 2,105

GO_annotations <- AnnotationDbi::select(GO.db,
                                        keys=unique(go_ipr$GO),
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>%
  dplyr::rename(TERM = GOID, TERM_NAME = TERM)

merged_ont <- go_ipr %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) %>%
  dplyr::filter(!is.na(TERM_NAME))

# BP
enGO_HDR_merged_BP <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_BP)

test <- as.data.table(enGO_HDR_merged_BP@result)

dotplot(enGO_HDR_merged_BP, showCategory = 40, title = "BP HDRs")

enGO_HDR_merged_plot_BP <- as.data.table(enGO_HDR_merged_BP@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(!is.na(Description)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(enGO_HDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6, color='black', hjust = 1),
        legend.text = element_text(size=5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.78, 0.42),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")

enGO_HDR_merged_plot <- as.data.table(enGO_HDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number()) %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen","oxidoreductase activity (2)", Description)) %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", "oxidoreductase activity (1)", Description)) 

plot_GO_MF <- ggplot(enGO_HDR_merged_plot) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.3, 3), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6, color='black', hjust = 1),
        legend.text = element_text(size=5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.74, 0.3),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        text = element_text(family="Helvetica"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(r = 10, b = 2, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF


final_plot <- cowplot::plot_grid(
  plot_ipr, plot_GO_BP, plot_GO_MF,
  rel_heights = c(1.5, 0.75, 0.9),
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("a","b","c"),
  label_size = 14,
  label_fontface = "bold")
final_plot

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/final_plotsAndData/IPR_domain_GO_BpMf_HDR_20251210.png", final_plot,  width = 7.5, height = 7, dpi = 600)












































#==============================================================================================================================================================================================================================#

# INTERPROSCAN (Gene Ontology) - just ARM genes as background - nonHDR arm domain genes

#==============================================================================================================================================================================================================================#
# BP
enGO_nHDR_merged_BP <- clusterProfiler::enricher(
  gene = nHD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_nHDR_merged_BP)

dotplot(enGO_nHDR_merged_BP, showCategory = 40, title = "BP nHDRs")

enGO_nHDR_merged_plot_BP <- as.data.table(enGO_nHDR_merged_BP@result) %>%
  tidyr::separate(GeneRatio, into = c("nhdr_gene_term", "nhdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (nhdr_gene_term / nhdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP_nHDR <- ggplot(enGO_nHDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_nHDR_merged_plot_BP$plotpoint, labels = enGO_nHDR_merged_plot_BP$Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_nHDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_nHDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_nHDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_nHDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_nHDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=9, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.35),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 15, r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP_nHDR



# MF
enGO_nHDR_merged_MF <- clusterProfiler::enricher(
  gene = nHD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_nHDR_merged_MF)

dotplot(enGO_nHDR_merged_MF, showCategory = 40, title = "MF nHDRs")

enGO_nHDR_merged_plot <- as.data.table(enGO_nHDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("nhdr_gene_term", "nhdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (nhdr_gene_term / nhdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_MF_nHDR <- ggplot(enGO_nHDR_merged_plot) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_nHDR_merged_plot$plotpoint, labels = enGO_nHDR_merged_plot$Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_nHDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_nHDR_merged_plot$Count, na.rm = TRUE) + min(enGO_nHDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_nHDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_nHDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=9, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.19),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(t = 15, r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF_nHDR


final_plot_nHDR <- cowplot::plot_grid(
  plot_GO_BP_nHDR, plot_GO_MF_nHDR,
  rel_heights = c(1.5, 3),
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("a", "b"),
  label_size = 14,
  label_fontface = "bold")
final_plot_nHDR

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/final_plotsAndData/IPR_GO_nHDR_armBackground.png", final_plot,  width = 7.5, height = 7.5, dpi = 600)
























































#==============================================================================================================================================================================================================================#

# INTERPROSCAN (Gene Ontology) - all genes as background - all HDR genes

#==============================================================================================================================================================================================================================#
IPR_GO_bckgrd <- unique(go_ipr$QX1410) # 11,756 genes

GO_annotations <- AnnotationDbi::select(GO.db,
                                        keys=unique(go_ipr$GO),
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>%
  dplyr::rename(TERM = GOID, TERM_NAME = TERM)

merged_ont <- go_ipr %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) %>%
  dplyr::filter(!is.na(TERM_NAME))

# BP
enGO_HDR_merged_BP <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_BP)

test <- as.data.table(enGO_HDR_merged_BP@result)

dotplot(enGO_HDR_merged_BP, showCategory = 40, title = "BP HDRs")

enGO_HDR_merged_plot_BP <- as.data.table(enGO_HDR_merged_BP@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(!is.na(Description)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(enGO_HDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=9, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.76, 0.4),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")

enGO_HDR_merged_plot <- as.data.table(enGO_HDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number()) %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (1)", Description))

plot_GO_MF <- ggplot(enGO_HDR_merged_plot) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=9, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.3),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF


#==============================================================================================================================================================================================================================#

# INTERPROSCAN (Gene Ontology) - all genes as background - nonHDR arm domain genes

#==============================================================================================================================================================================================================================#
go_ipr_arms <- go_ipr %>% dplyr::filter(QX1410 %in% arm_genes)
IPR_GO_bckgrd_arms <- unique(go_ipr_arms$QX1410) # 3,900 genes

GO_annotations <- AnnotationDbi::select(GO.db,
                                        keys=unique(go_ipr$GO),
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>%
  dplyr::rename(TERM = GOID, TERM_NAME = TERM)

merged_ont <- go_ipr %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) %>%
  dplyr::filter(!is.na(TERM_NAME))

# BP
enGO_HDR_merged_BP <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_BP)

test <- as.data.table(enGO_HDR_merged_BP@result)

dotplot(enGO_HDR_merged_BP, showCategory = 40, title = "BP HDRs")

enGO_HDR_merged_plot_BP <- as.data.table(enGO_HDR_merged_BP@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(!is.na(Description)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number())

plot_GO_BP <- ggplot(enGO_HDR_merged_plot_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot_BP$plotpoint, labels = enGO_HDR_merged_plot_BP$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot_BP$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=9, color='black', face = 'bold'),
        axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.76, 0.4),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(b = 5, r = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP



# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_annotations %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")

enGO_HDR_merged_plot <- as.data.table(enGO_HDR_merged_MF@result) %>%
  tidyr::separate(GeneRatio, into = c("hdr_gene_term", "hdr_gene_hit"), sep = "/", convert = TRUE) %>%
  tidyr::separate(BgRatio, into = c("bgd_term", "bgd_hit"), sep = "/", convert = TRUE) %>%
  dplyr::mutate(EnrichRatio = (hdr_gene_term / hdr_gene_hit) / (bgd_term / bgd_hit)) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::mutate(plotpoint = dplyr::row_number()) %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (1)", Description))

plot_GO_MF <- ggplot(enGO_HDR_merged_plot) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = enGO_HDR_merged_plot$plotpoint, labels = enGO_HDR_merged_plot$Description, name = "", expand = c(0.06,0.06)) +
  scale_fill_gradient(low = "blue", high = "red", breaks = c(round(min(enGO_HDR_merged_plot$Count, na.rm = TRUE)), round((max(enGO_HDR_merged_plot$Count, na.rm = TRUE) + min(enGO_HDR_merged_plot$Count, na.rm = TRUE) ) / 2), round(max(enGO_HDR_merged_plot$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(enGO_HDR_merged_plot$EnrichRatio, n = 3)) +
  theme(axis.text.x = element_text(size=8, color='black'),
        axis.text.y = element_text(size=8, color='black'),
        axis.title = element_text(size=9, color='black', face = 'bold'),
        # plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_blank(),
        legend.title = element_text(size=6.5, color='black', hjust = 1),
        legend.text = element_text(size=5.5, color='black', hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.3),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.0001, 'cm'),
        legend.key.height = unit(0.01, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.box.just = "right",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.margin = margin(r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF








#==============================================================================================================================================================================================================================#

# Look at Orthology of N2 genes previoulsy shown to be enriched in environmental and pathogenic responses

#==============================================================================================================================================================================================================================#
# Looking at monooxygenase GO MF term enriched in HDR arm genes with arm genes background
monoOx <- enGO_HDR_merged_plot %>% dplyr::filter(Description == "monooxygenase activity") %>% dplyr::select(geneID) %>% tidyr::separate_rows(geneID, sep = '/') %>% dplyr::pull()

monoOxOrthos <- ortho_genes_dd %>%
  dplyr::select(Orthogroup, QX1410, N2) %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::filter(QX1410 %in% monoOx) 

monoOx_N2 <- monoOxOrthos %>%
  dplyr::select(N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ", ") 

N2_name <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/misc/N2_WBGeneID_Name.tsv", col_names = F) %>% dplyr::rename(WBGeneID = X1, name = X2)

N2_ortho_name <- monoOx_N2 %>%
  dplyr::left_join(N2_name, by = c("N2" = "WBGeneID"))

# Looking at the innate immune response GO BP term enriched in HDR arm genes with arm genes background

monoIIR <- enGO_HDR_merged_plot_BP %>% dplyr::filter(Description == 'innate immune response') %>% dplyr::select(geneID) %>% tidyr::separate_rows(geneID, sep = '/') %>% dplyr::pull()

monoIIR_orthos <- ortho_genes_dd %>%
  dplyr::select(Orthogroup, QX1410, N2) %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::filter(QX1410 %in% monoIIR)

monoIIR_N2 <- monoIIR_orthos %>%
  dplyr::select(N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ', ')

N2_IIR_ortho_names <- monoIIR_N2 %>%
  dplyr::left_join(N2_name, by = c("N2" = "WBGeneID"))



# Looking at xenobiotic metabolic process
xenoIPR <- enGO_HDR_merged_plot_BP %>% dplyr::filter(Description == 'xenobiotic metabolic process') %>% dplyr::select(geneID) %>% tidyr::separate_rows(geneID, sep = '/') %>% dplyr::pull()

xenoIPR_orthos <- ortho_genes_dd %>%
  dplyr::select(Orthogroup, QX1410, N2) %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::filter(QX1410 %in% xenoIPR)

xenoIPR_N2 <- xenoIPR_orthos %>%
  dplyr::select(N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ', ')

N2_IPR_ortho_names <- xenoIPR_N2 %>%
  dplyr::left_join(N2_name, by = c("N2" = "WBGeneID"))

