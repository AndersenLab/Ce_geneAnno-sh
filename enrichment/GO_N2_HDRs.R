#library(tidyverse)
library(plyr)
library(valr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(data.table)
library(org.Ce.eg.db)
# BiocManager::install("org.Ce.eg.db")
#library(biomaRt) ##BiocManager::install("biomaRt")
# install.packages("ggfun")
# install.packages("remotes")
# remotes::install_github("YuLab-SMU/ggtree")
library(clusterProfiler) ## BiocManager::install("clusterProfiler")
library(enrichplot)
# BiocManager::install("GenomicRanges")
library(GenomicRanges)
# BiocManager::install("rtracklayer")
library(rtracklayer)


### GO enrichment tests
enrich_go <- function(x = NULL){
  
  # GO term enrichment  
  wb_ids <- dplyr::filter(x) %>% dplyr::pull(WBGeneID)
  
  mf <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "bonferroni",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  bp <- enrichGO(gene          = wb_ids,
                 OrgDb         = org.Ce.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "bonferroni",
                 keyType       = 'ENSEMBL',
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
  
  return(list(mf, bp))
}


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

gffCE <-  ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/misc/N2.WBonly.WS283.PConly.gff3")

gffFinal <- gffCE %>%
  dplyr::filter(type=="gene") %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes, into=c("ID","rest"), sep = ";Name=") %>%
  dplyr::mutate(ID=gsub("ID=Gene:","",ID)) %>%
  dplyr::select(seqid,start,end,ID) %>%
  dplyr::rename(CHROM=seqid)

nHDR_armDomain <- readr::read_csv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/elegans/GO_terms/chromosome_domain_Celegans.csv") %>%
  dplyr::rename(Chromosome=chrom) %>%
  dplyr::mutate(left=left*1e3,right=right*1e3) 

hdr_regions <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("CHROM","minStart","maxEnd","STRAIN")) 
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



HDreg_arm <- ldply(reglist,data.frame)
HDreg_center <- ldply(reglist_center,data.frame)
HDgene <- getGenesFromBed(gffFinal,HDreg_arm)
HD_gene_vector <- HDgene[[]]
HD_N2_genes <- (do.call(rbind, HDgene))
HD_gene_vector <- unique(HD_N2_genes$N2) 

nHDreg <- classifyingArms(HDreg_arm %>% dplyr::rename(minStart=start,maxEnd=end), nHDR_armDomain) %>% dplyr::select(-domain)
nHDgene <- getGenesFromBed(gffFinal,nHDreg)
nHD_N2_genes <- do.call(rbind, nHDgene)
nHD_gene_vector <- unique(nHD_N2_genes$N2) 


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
  geom_rect(data=HD_N2_genes, aes(xmin=start,xmax=end,ymin=-1,ymax=1)) +
  facet_wrap(~CHROM,ncol=1,scales="free_x") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  xlab("Physical postion (Mb)")
posplot2


### Go enrichment test for the divergent regions
gene_files_arm <- ldply(nHDgene,data.frame) %>%
  dplyr::select(ID) %>%
  dplyr::rename(WBGeneID=ID)

df_GO_nondiv_arm <- enrich_go(as.data.table(gene_files_arm))

gene_name_Pro_BP_arm <- setReadable(df_GO_nondiv_arm[[2]], OrgDb = org.Ce.eg.db)
# barplot(gene_name_Pro_BP_arm, showCategory=100)

gene_name_Pro_MF_arm <- setReadable(df_GO_nondiv_arm[[1]], OrgDb = org.Ce.eg.db)
# barplot(gene_name_Pro_MF_arm, showCategory=100)

## all 
gene_files_div_all <- ldply(HDgene,data.frame) %>% 
  dplyr::select(ID) %>%
  dplyr::rename(WBGeneID=ID)

df_GO_div_all <- enrich_go(as.data.table(gene_files_div_all))

gene_name_Pro_BP_div_all <- setReadable(df_GO_div_all[[2]], OrgDb = org.Ce.eg.db)
#barplot(gene_name_Pro_BP_div_all, showCategory=100) 

gene_name_Pro_MF_div_all <- setReadable(df_GO_div_all[[1]], OrgDb = org.Ce.eg.db)
#barplot(gene_name_Pro_MF_div_all)



### Biological process (BP)
df_GO_enrich_BP <- rbind(dplyr::mutate(gene_name_Pro_BP_arm[], freq="control_arms"), dplyr::mutate(gene_name_Pro_BP_div_all[], freq="all")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(p.adjust)

GO_list_BP <- df_GO_enrich_BP$Description

GO_list_BP_plotpoint <- data.frame(Description=GO_list_BP, plotpoint=length(GO_list_BP):1)

df_GO_enrich_BP_sum <- df_GO_enrich_BP %>%
  dplyr::filter(Description %in% GO_list_BP) %>%
  dplyr::left_join(., GO_list_BP_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_BP <- df_GO_enrich_BP_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

df_GO_enrich_BP_sum$freq <- factor(df_GO_enrich_BP_sum$freq, levels = c("all","control_arms"), labels = c("Hyper-divergent regions","Non-divergent arms"))

plot_GO_BP <- df_GO_enrich_BP_sum %>%
  ggplot(.) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x=-log10(p.adjust), y=plotpoint, size=EnrichRatio, fill=Count), shape=21) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_BP):1, labels = GO_list_BP, name = "",sec.axis = dup_axis(labels = df_class_total_BP$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "tan1", high = "#DB6333",breaks = c(50, 100, 150)) +
  scale_shape_manual(values=c(21, 22)) +
  scale_size_continuous(range = c(2,10), breaks = c(2.5,3,4)) +
  theme(axis.text = element_text(size=12, color='black'), 
        axis.title = element_text(size=12, color='black'), 
        legend.title = element_text(size=12, color='black'), 
        legend.text = element_text(size=12, color='black'), 
        #legend.position = "inside", 
        legend.position = "inside", 
        legend.position.inside = c(0.78, 0.12),
        panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(
         # shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
         size= guide_legend(nrow=1, order = 2, title.position = "top"), 
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="-log10 (corrected p-value)", size = "Fold enrichment", fill = "Gene counts") # shape = "Frequency"


plot_GO_BP


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/elegans/GOBP_HDRs.png", plot_GO_BP,height = 8.5,width = 12.5,units = 'in',device = 'png',dpi = 900)


# new_term <- df_GO_enrich_BP_sum %>%
#   dplyr::filter(Description == "biological process involved in interspecies interaction between organisms") %>%
#   dplyr::select(ID,Description,geneID,Count) %>%
#   dplyr::mutate(geneID = gsub("/",",",geneID))
# 
# print(new_term$geneID)



### Molecular function (MF)
# 
df_GO_enrich_MF <- rbind(dplyr::mutate(gene_name_Pro_MF_arm[], freq="control_arms"), dplyr::mutate(gene_name_Pro_MF_div_all[], freq="all")) %>%
  tidyr::separate(GeneRatio, c("N_div_class", "N_div_all")) %>%
  tidyr::separate(BgRatio, c("Total_class", "Total_genome")) %>%
  dplyr::mutate(N_div_class = as.numeric(N_div_class), N_div_all = as.numeric(N_div_all), Total_class = as.numeric(Total_class), Total_genome = as.numeric(Total_genome)) %>%
  dplyr::mutate(GeneRatio = N_div_class/N_div_all, BgRatio = Total_class/Total_genome) %>%
  dplyr::mutate(EnrichRatio = GeneRatio/BgRatio) %>%
  dplyr::arrange(p.adjust)

df_GO_enrich_MF$Description <- gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen", "oxidoreductase activity", df_GO_enrich_MF$Description)
df_GO_enrich_MF$Description <- gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen", "oxidoreductase activity (b)", df_GO_enrich_MF$Description)

GO_list_MF <- df_GO_enrich_MF$Description

GO_list_MF_plotpoint <- data.frame(Description=GO_list_MF, plotpoint=length(GO_list_MF):1)

df_GO_enrich_MF_sum <- df_GO_enrich_MF %>%
  dplyr::filter(Description %in% GO_list_MF) %>%
  dplyr::left_join(., GO_list_MF_plotpoint, by = "Description") %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(class_gene_total = max(Total_class)) %>%
  dplyr::ungroup()

df_class_total_MF <- df_GO_enrich_MF_sum %>%
  dplyr::distinct(Description, class_gene_total, plotpoint) %>%
  dplyr::arrange(-plotpoint)

df_GO_enrich_MF_sum$freq <- factor(df_GO_enrich_MF_sum$freq, levels = c("all","control_arms"), labels = c("Divergent regions","Non-divergent arms"))

plot_GO_MF <- df_GO_enrich_MF_sum %>%
  ggplot(.) +
  geom_vline(xintercept = -log10(0.05), color='blue', size=0.4) +
  geom_point(alpha=0.9) +
  aes(x=-log10(p.adjust), y=plotpoint, shape="21", size=EnrichRatio, fill=Count, alpha=0.9) +
  theme_bw() +
  #scale_x_continuous(breaks = c(0,2,4,6,8,10), expand = c(0.2,0.2)) +
  scale_y_continuous(breaks = length(GO_list_MF):1, labels = GO_list_MF, name = "", expand = c(0.05,0.05)) +
  #scale_y_continuous(breaks = length(GO_list_MF):1, labels = GO_list_MF, name = "",sec.axis = dup_axis(labels = df_class_total_MF$class_gene_total, name = "Total gene counts")) +
  scale_fill_gradient(low = "lightpink", high = "red") +
  scale_shape_manual(values=21) +
  scale_size_continuous(range = c(2,4), breaks = c(2,3,4)) +
  theme(axis.text = element_text(size=7, color='black'),
        axis.title = element_text(size=7, color='black'),
        legend.title = element_text(size=7, color='black'),
        legend.text = element_text(size=7, color='black'),
        legend.position = c(-1,0.05),
        panel.grid = element_blank(),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(size= guide_legend(nrow=1, order = 2, title.position = "top"),
         fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 2)) +
  labs(x="-log10(Corrected P-value)", shape = "", size = "Fold enrichment", fill = "Gene counts")

plot_GO_MF








################## BOTH BP AND MF #######################
bp_mf <- df_GO_enrich_MF_sum %>% dplyr::bind_rows(df_GO_enrich_BP_sum) %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 20) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::mutate(plotpoint = row_number())

hello <- ggplot(bp_mf) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x=-log10(p.adjust), y=plotpoint, size=EnrichRatio, fill=Count), shape=21) +
  theme_bw() +
  scale_fill_gradient(low = "tan1", high = "#DB6333",breaks = c(50, 100, 150)) +
  scale_shape_manual(values=c(21, 22)) +
  scale_y_continuous(
    breaks = bp_mf$plotpoint,
    labels = bp_mf$Description  # or whatever your GO term column is called
  ) +
  scale_size_continuous(range = c(4,15), breaks = c(2.5,3,4)) +
  theme(axis.text = element_text(size=16, color='black'), 
        axis.title = element_text(size=18, color='black'), 
        axis.title.y = element_blank(),
        legend.title = element_text(size=20, color='black'), 
        legend.text = element_text(size=20, color='black'), 
        #legend.position = "inside", 
        legend.position = "inside", 
        legend.position.inside = c(0.78, 0.12),
        panel.grid = element_blank(),
        legend.direction = "horizontal", legend.box = "vertical",
        legend.spacing.y = unit(0.05, 'in')) +
  guides(
    # shape= guide_legend(nrow=1, order = 3, title.position = "top", override.aes = list(size=3)), 
    size= guide_legend(nrow=1, order = 2, title.position = "top"), 
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", label.vjust = 1)) +
  labs(x="-log10 (corrected p-value)", size = "Fold enrichment", fill = "Gene counts") # shape = "Frequency"
hello





# 
# plot_GO_BP_MF <- cowplot::plot_grid(plot_GO_BP, plot_GO_MF, ncol=1, rel_heights = c(1,1), labels = c("A", "B"), align = "HD", label_size = 12)
# 
# ggsave("/projects/b1059/projects/Nicolas/hyperdivergent_regions/elegans/final_f4_GO.png",plot_GO_BP_MF,height = 5.5,width = 6.5,units = 'in',device = 'png',dpi = 900) 
# 
# #ggsave("Plots/supplementary_figures/plot_GO_BP_MF.pdf", plot = plot_GO_BP_MF, width = 7.5, height =9, units = "in")
# 
# #save(GO_list_BP, GO_list_BP_plotpoint, df_GO_enrich_BP_sum, GO_list_MF, GO_list_MF_plotpoint, df_GO_enrich_MF_sum, file="Processed_Data/GO_enrich.RData")
