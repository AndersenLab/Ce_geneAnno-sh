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
#library(broom) # aded by Nic


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
# GO - extract any-to-ones between QX and N2 and lift over GO terms #
# ======================================================================================================================================================================================== #
# Extracting QX - N2 any-to-ones
qn <- all_relations %>%
  dplyr::select(-AF16_count) %>%
  dplyr::filter(QX1410_count == 1) %>%  ################################# FILTERING FOR ONE-TO-ONES 
  # dplyr::filter(QX1410_count >= 1) %>%
  dplyr::filter(N2_count == 1) #%>% 
  # dplyr::filter(HOG != "N0.HOG0000062") %>%
  # dplyr::filter(HOG != "N0.HOG0000172")

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

# Do not use GO IDs from N2 htis... use eggNOG on remaining QX1410 genes that were not able to have GO IDs lifted over from N2 HOGs, N2 RBBH, and AF16 HOGs 
# N2_euk_bgd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/eggNOG_mapper/elegans/output/euk_annotations/N2_longestIso_background_euk_GOannotations.tsv", col_names = c("transcript","GO", "term_name")) %>%
# dplyr::left_join(tran_gene, by = "transcript") %>%
# tidyr::separate_rows(GO, sep = ",") 

n2_gene_database <- gene2go %>%
  dplyr::select(WORMBASE, GO, ONTOLOGY) %>%
  dplyr::mutate(n2_gene = ifelse(WORMBASE %in% n2_wbGenes$gene, TRUE, FALSE)) %>%
  dplyr::mutate(n2_gene = ifelse(n2_gene == TRUE, WORMBASE, n2_gene)) 

merged  <- n2_gene_database %>%
  dplyr::select(n2_gene,GO) %>%
  dplyr::filter(n2_gene != FALSE) %>% 
  dplyr::filter(!is.na(GO)) %>%
  dplyr::rename(gene = n2_gene) # 14,033 N2 genes


# Lift over GO terms from N2 to QX1410
GO_brig <- oneoneqn %>%
  tidyr::separate_rows(QX1410, sep = ",\\s*") %>% 
  dplyr::left_join(merged, by = c("N2" = "gene"), relationship = "many-to-many") %>%
  dplyr::filter(!is.na(GO)) # 10,749


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
  dplyr::left_join(rbbh, by = "N2") #%>%
  # dplyr::filter(HOG != "N0.HOG0000062") %>% # gene model error in this HOG - MANY QX1410 genes for a single N2 gene because of some AA similarity for different overlapping genes on Watson and Crick
  # dplyr::filter(HOG != "N0.HOG0000172")

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
  dplyr::filter(!is.na(GO_final)) # 11,334 QX1410 genes



# ======================================================================================================================================================================================== #
# After adding output from RBBH, then look at QX1410 - AF16 any-to-ones and then convert AF16 gene name to N2 WBGeneID and pull any GO IDs #
# ======================================================================================================================================================================================== #
cb_alias <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/CB_N2_alias.tsv", col_names = c("AF_WBGeneID","alias"))
N2_alias <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/N2_aliasNames.tsv", col_names = c("N2_WBGeneID","alias")) %>%
  dplyr::select(alias,N2_WBGeneID)


# Isolating QX1410 - AF16 any-to-one
cbcb <- all_relations %>%
  dplyr::select(-N2_count) %>%
  dplyr::filter(QX1410_count == 1) %>% ################################# FILTERING FOR ONE-TO-ONES 
  # dplyr::filter(QX1410_count >= 1) %>%
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
  dplyr::filter(HOG %in% cbcb$HOG) %>% ################################# FILTERING FOR ONE-TO-ONES 
  dplyr::select(QX1410) %>%
  tidyr::separate_rows(QX1410, sep = ',\\s') %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::left_join(join_final, by = "QX1410") %>%
  dplyr::left_join(cbcb_genes, by = "QX1410") %>%
  dplyr::select(QX1410, N2.x, GO_final.x, N2_WBGeneID, GO_fromAlias) %>%
  dplyr::mutate(N2 = ifelse(!is.na(N2.x), N2.x, N2_WBGeneID)) %>%
  dplyr::mutate(GO = ifelse(!is.na(GO_final.x), GO_final.x, GO_fromAlias)) %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::select(QX1410, GO) # 11,568 

# 11,363 one-to-one's
# validating.....
oneone_HOG_list <- ortho_genes_dd %>%
  dplyr::filter(HOG %in% cbcb$HOG) %>% ################################# FILTERING FOR ONE-TO-ONES 
  dplyr::select(HOG,QX1410) %>%
  tidyr::separate_rows(QX1410, sep = ',\\s') %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::left_join(join_final, by = "QX1410") %>%
  dplyr::left_join(cbcb_genes, by = "QX1410") %>%
  dplyr::select(HOG,QX1410, N2.x, GO_final.x, N2_WBGeneID, GO_fromAlias) %>%
  dplyr::mutate(N2 = ifelse(!is.na(N2.x), N2.x, N2_WBGeneID)) %>%
  dplyr::mutate(GO = ifelse(!is.na(GO_final.x), GO_final.x, GO_fromAlias)) %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::select(HOG) %>%
  dplyr::distinct() %>%
  dplyr::pull()

filtered <- all_relations %>%
  dplyr::filter(HOG %in% oneone_HOG_list)




# Quantifying how many one-to-many N2-QX1410 relationships there are
count_complexHOGs <- ortho_genes_dd %>%
  dplyr::select(HOG,OG,QX1410,N2,AF16) %>%
  dplyr::filter(!is.na(N2), !is.na(QX1410)) #%>%
  # dplyr::filter(HOG != "N0.HOG0000062", HOG != "N0.HOG0000172") # already validated incorrect HOGs because of gene modeling errors resulting in fusions

gene_count <- count_complexHOGs
col_names <- colnames(count_complexHOGs)

for (i in 3:length(col_names)) {
  print(paste0(i,"out of", length(col_names)))
  temp_colname = paste0(col_names[i], "_count")
  
  gene_count <- gene_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(col_names[i]),", ") + 1)
}

howMany <- gene_count %>%
  dplyr::mutate(complex = ifelse(N2_count == 1 & QX1410_count >=2, T,F)) %>%
  dplyr::mutate(ratio = ifelse(complex == T, QX1410_count / N2_count, NA))

ggplot(howMany) +
  geom_histogram(aes(x = ratio), binwidth = 1, fill = 'purple') +
  theme(
    axis.title = element_text(size = 14, face = 'bold', color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = NA)) +
  scale_x_continuous(expand = c(0.01,0)) +
  scale_y_continuous(expand = c(0.01,0)) +
  xlab("N2 genes / QX1410 genes") +
  ylab("Count")


# Assessing liftover for WormCat terms and N2 genes associated with them for Nic - 07/11 #
# test <- ortho_genes_dd %>%
#   dplyr::select(QX1410) %>%
#   tidyr::separate_rows(QX1410, sep = ',\\s') %>%
#   dplyr::filter(!is.na(QX1410)) %>%
#   dplyr::left_join(join_final, by = "QX1410") %>%
#   dplyr::left_join(cbcb_genes, by = "QX1410") %>%
#   dplyr::select(QX1410, N2.x, GO_final.x, N2_WBGeneID, GO_fromAlias) %>%
#   dplyr::mutate(N2 = ifelse(!is.na(N2.x), N2.x, N2_WBGeneID)) %>%
#   dplyr::mutate(GO = ifelse(!is.na(GO_final.x), GO_final.x, GO_fromAlias)) %>%
#   dplyr::filter(!is.na(GO)) %>%
#   dplyr::select(QX1410,GO,N2) 
# 
# pulled <- test %>%
#   dplyr::select(N2) %>%
#   dplyr::distinct() %>%
#   dplyr::pull()

# terms <- readr::read_csv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/misc/N2_WormCat_terms_andGenes.csv", col_names = c("seqID", "N2", "class","cat1","cat2", "cat3", "desc")) 
# 
# gpcr <- terms %>%
#   dplyr::filter(grepl("Signaling: heteromeric G protein", cat1)) %>%
#   dplyr::select(N2, cat1) %>%
#   dplyr::mutate(lifted_over = ifelse(N2 %in% pulled, TRUE, FALSE)) 
# summary_gp <- gpcr %>%
#   dplyr::count(lifted_over)
# test_gpcr <- test %>%
#   dplyr::left_join(gpcr, by = "N2") %>%
#   dplyr::mutate(inHDR = ifelse(QX1410 %in% HD_all_gene_vector, TRUE, FALSE)) %>%
#   dplyr::filter(!is.na(lifted_over)) %>%
#   dplyr::distinct(QX1410, .keep_all = T) %>%
#   dplyr::count(inHDR)
# summary_gpHDR <- test_gpcr %>%
#   dplyr::count(inHDR)
# 
#   
# Clectin <- terms %>%
#   dplyr::filter(grepl("Stress response: C-type Lectin", cat1)) %>%
#   dplyr::select(N2, cat1) %>%
#   dplyr::mutate(lifted_over = ifelse(N2 %in% pulled, TRUE, FALSE))
# summary_Cl <- Clectin %>%
#   dplyr::count(lifted_over)
# test_Clectin <- test %>%
#   dplyr::left_join(Clectin, by = "N2") %>%
#   dplyr::mutate(inHDR = ifelse(QX1410 %in% HD_all_gene_vector, TRUE, FALSE)) %>%
#   dplyr::filter(!is.na(lifted_over)) %>%
#   dplyr::distinct(QX1410, .keep_all = T) %>%
#   dplyr::count(inHDR)
#                 
# 
# E3 <- terms %>%
#   dplyr::filter(grepl("Proteolysis proteasome: E3", cat1)) %>%
#   dplyr::select(N2, cat1) %>%
#   dplyr::mutate(lifted_over = ifelse(N2 %in% pulled, TRUE, FALSE)) 
# summary_E3 <- E3 %>%
#   dplyr::count(lifted_over)
# test_E3 <- test %>%
#   dplyr::left_join(E3, by = "N2") %>%
#   dplyr::mutate(inHDR = ifelse(QX1410 %in% HD_all_gene_vector, TRUE, FALSE)) %>%
#   dplyr::filter(!is.na(lifted_over)) %>%
#   dplyr::distinct(QX1410, .keep_all = T) %>%
#   dplyr::count(inHDR)
# 
#                
# Fbox <- terms %>% 
#   dplyr::filter(grepl("Proteolysis proteasome: E3: F box", cat2)) %>%
#   dplyr::select(N2, cat2) %>%
#   dplyr::mutate(lifted_over = ifelse(N2 %in% pulled, TRUE, FALSE)) 
# summary_Fbox <- Fbox %>%
#   dplyr::count(lifted_over)
# test_Fbox <- test %>%
#   dplyr::left_join(Fbox, by = "N2") %>%
#   dplyr::mutate(inHDR = ifelse(QX1410 %in% HD_all_gene_vector, TRUE, FALSE)) %>%
#   dplyr::filter(!is.na(lifted_over)) %>%
#   dplyr::distinct(QX1410, .keep_all = T) %>%
#   dplyr::count(inHDR)
# 
# pals <- N2_alias %>%
#   dplyr::filter(grepl("pals", alias)) %>%
#   dplyr::rename(N2 = N2_WBGeneID) %>%
#   dplyr::mutate(lifted_over = ifelse(N2 %in% pulled, TRUE, FALSE))
# test_pals <- test %>%
#   dplyr::left_join(pals, by = "N2") %>%
#   dplyr::mutate(inHDR = ifelse(QX1410 %in% HD_all_gene_vector, TRUE, FALSE)) %>%
#   dplyr::filter(!is.na(lifted_over)) %>%
#   dplyr::distinct(QX1410, .keep_all = T) #%>%
#   dplyr::count(inHDR)



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

hdr_regions <- readr::read_tsv("/vast/eande106/projects/Nicolas/hyperdivergent_regions/briggsae/multi_reference/reference_realign/HDR_CB_allStrain_5kbclust_20250730.tsv") %>%
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


HDreg_arm <- ldply(reglist,data.frame)
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

nHDreg <- classifyingArms(HDreg_arm %>% dplyr::rename(minStart=start,maxEnd=end), nHDR_armDomain) %>% dplyr::select(-domain)
nHDgene <- getGenesFromBed(gffFinal,nHDreg)
nHD_QX_genes <- do.call(rbind, nHDgene)
nHD_gene_vector <- unique(nHD_QX_genes$QX1410) 


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


arm_lifted <- all_arm_genes_df %>% dplyr::filter(QX1410 %in% final_GO$QX1410)
all_lifted <- gffFinal %>% dplyr::filter(QX1410 %in% final_GO$QX1410)

bloop <- ggplot(data = nHDR_armDomain) + 
  geom_rect(aes(xmin = left / 1e6, xmax = right / 1e6, ymin = 0.5, ymax = 1.5), fill = 'black') +
  geom_rect(data = arm_lifted %>% dplyr::rename(Chromosome = CHROM), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.5, ymax = 1.5), fill = 'gold') +
  geom_rect(data = all_lifted %>% dplyr::rename(Chromosome = CHROM), aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0.75, ymax = 1.25), fill = "skyblue") +
  facet_wrap(~Chromosome, scales = 'free') + 
  theme(
    legend.position = 'none')
bloop 

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/misc/lift_overPlot.png", bloop, width = 14, height= 8, dpi = 300)

QX_arm_genes <- final_GO %>%
  dplyr::filter(QX1410 %in% arm_genes)
# QX1410_merged_bgd_genes <- unique(final_GO$QX1410) #  12,532 genes - 9,415 when filtering for one-to-ones
QX1410_merged_bgd_genes <- unique(QX_arm_genes$QX1410) #  12,532 genes - 9,415 when filtering for one-to-ones - 2,518 when restricted to genes only found in arm domains

# Dumb test to ensure that it doesn't matter if we use ONLY genes in HDR arms that have GO ID annotation (shouldn't matter because the enrichmenet ratio is dependent on the ones that are found in the DB)
# blah <- HD_QX_genes %>% dplyr::select(QX1410) %>% dplyr::filter(QX1410 %in% QX1410_merged_bgd_genes)
# HD_gene_vector <- dplyr::pull(blah)
GO_terms_merged <- AnnotationDbi::select(GO.db, 
                                         keys=unique(final_GO$GO),   
                                         columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                         keytype="GOID") %>% 
  dplyr::rename(TERM = GOID, TERM_NAME = TERM) # 8,671

merged_ont <- final_GO %>%
  dplyr::left_join(GO_terms_merged, by = c("GO" = "TERM")) 









############################################################
# All arm genes against all QX1410 genes that received GO ID liftover
# QX_arm_genes <- final_GO %>%
#   dplyr::filter(QX1410 %in% arm_genes)
# QX1410_merged_bgd_genes <- unique(final_GO$QX1410) #  12,532 genes - 9,415 when filtering for one-to-ones
# arm_genes <- unique(QX_arm_genes$QX1410) #  12,532 genes - 9,415 when filtering for one-to-ones - 2,518 when restricted to genes only found in arm domains
# 
# # Dumb test to ensure that it doesn't matter if we use ONLY genes in HDR arms that have GO ID annotation (shouldn't matter because the enrichmenet ratio is dependent on the ones that are found in the DB)
# # blah <- HD_QX_genes %>% dplyr::select(QX1410) %>% dplyr::filter(QX1410 %in% QX1410_merged_bgd_genes)
# # HD_gene_vector <- dplyr::pull(blah)
# GO_terms_merged <- AnnotationDbi::select(GO.db, 
#                                          keys=unique(final_GO$GO),   
#                                          columns = c("TERM", "DEFINITION", "ONTOLOGY"),
#                                          keytype="GOID") %>% 
#   dplyr::rename(TERM = GOID, TERM_NAME = TERM) # 8,671
# 
# merged_ont <- final_GO %>%
#   dplyr::left_join(GO_terms_merged, by = c("GO" = "TERM"))
# 
# 
# enGO_HDR_merged <- clusterProfiler::enricher(
#   gene = arm_genes,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX1410_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_HDR_merged)
# 
# dotplot(enGO_HDR_merged, showCategory = 40, title = "BP arm genes")
# 
# # MF
# enGO_HDR_merged_MF <- clusterProfiler::enricher(
#   gene = arm_genes,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX1410_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_HDR_merged_MF)
# 
# dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF arm genes")











### Enrichment for HDR arm genes against arm genes background
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

GO_list_BP <- GO_list_BP %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::mutate(Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)), plotpoint = dplyr::row_number())

# Dimensions for figure #
plot_GO_BP <- ggplot(GO_list_BP) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_BP$plotpoint, labels = GO_list_BP$Wrapped_Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(GO_list_BP$Count, na.rm = TRUE)), round((max(GO_list_BP$Count, na.rm = TRUE) + min(GO_list_BP$Count, na.rm = TRUE) ) / 2), round(max(GO_list_BP$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5,4), name = "Fold enrichment", breaks = pretty(GO_list_BP$EnrichRatio, n = 3)) +
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
        plot.margin = margin(t = 10, r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:BP terms for genes in HDRs", x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_BP


GOBP_tableSave <- enGO_HDR_dt %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::select(-plotpoint) %>%
  dplyr::mutate(geneID = gsub("/", ", ", geneID))

# Need to revert-back any simplifed Descriptions
# write.table(GOBP_tableSave, file = "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOBP_HDRgenes.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Adding N2 genes, aliases, HOGs, and OGs
qn_table <- ortho_genes_dd %>%
  dplyr::select(QX1410) %>%
  tidyr::separate_rows(QX1410, sep = ',\\s') %>%
  dplyr::filter(!is.na(QX1410)) %>%
  dplyr::left_join(join_final, by = "QX1410") %>%
  dplyr::left_join(cbcb_genes, by = "QX1410") %>%
  dplyr::select(QX1410, N2.x, GO_final.x, N2_WBGeneID, GO_fromAlias) %>%
  dplyr::mutate(N2 = ifelse(!is.na(N2.x), N2.x, N2_WBGeneID)) %>%
  dplyr::mutate(GO = ifelse(!is.na(GO_final.x), GO_final.x, GO_fromAlias)) %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::select(QX1410, N2) %>%
  dplyr::distinct()


n2_alias <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/N2_alias.tsv", col_names = c("N2_WBGeneID","alias")) %>%
  dplyr::select(alias,N2_WBGeneID)

n2_hogs_ogs <- ortho_genes_dd %>%
  dplyr::select(HOG,OG,N2) %>%
  dplyr::filter(!is.na(N2)) %>%
  tidyr::separate_rows(N2, sep = ",\\s") %>%
  dplyr::select(N2,OG,HOG)

N2_tableSave <- GOBP_tableSave %>%
  tidyr::separate_rows(geneID, sep = ",\\s") %>%
  dplyr::left_join(qn_table, by = c("geneID" = "QX1410")) %>%
  dplyr::left_join(n2_alias, by = c("N2" = "N2_WBGeneID")) %>%
  dplyr::arrange(Description, N2) %>%
  dplyr::left_join(n2_hogs_ogs, by = "N2") %>%
  dplyr::group_by(ID,Description,hdr_gene_term,hdr_gene_hit,bgd_term,bgd_hit,pvalue,p.adjust,qvalue,Count,EnrichRatio) %>%
  dplyr::summarise(
    geneID = paste(unique(geneID), collapse = ","),
    N2 = paste(unique(N2), collapse = ","),
    alias = paste(unique(alias), collapse = ","),
    HOGs = paste(unique(HOG), collapse = ","),
    OGs = paste(unique(OG), collapse = ","),
    .groups = "drop")

# write.table(N2_tableSave, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/final_plotsAndData/GOBP_HDRgenes_updated_NOsulp7.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# MF
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
  # dplyr::slice_head(n=30) %>%
  dplyr::mutate(Wrapped_Description = insert_linebreaks(Description),
                Wrapped_Description = factor(Wrapped_Description, levels = unique(Wrapped_Description)),
                plotpoint = dplyr::row_number()) 

GO_list_MF <- GO_list_MF %>%
  dplyr::mutate(Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen","oxidoreductase activity (1)", Description),
                Description = gsub("oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen","oxidoreductase activity (2)", Description),
                Description = gsub("transmitter-gated ion channel activity involved in regulation of postsynaptic membrane potential", "transmitter-gated ion channel (1)", Description))

# ### Dimensions for figure ###
plot_GO_MF <- ggplot(GO_list_MF) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(p.adjust), y = plotpoint, size = EnrichRatio, fill = Count), shape = 21) +
  scale_y_continuous(breaks = GO_list_MF$plotpoint, labels = GO_list_MF$Description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(GO_list_MF$Count, na.rm = TRUE)), round((max(GO_list_MF$Count, na.rm = TRUE) + min(GO_list_MF$Count, na.rm = TRUE) ) / 2), round(max(GO_list_MF$Count, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(GO_list_MF$EnrichRatio, n = 3)) +
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
        plot.margin = margin(t = 10, r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched GO:MF terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_GO_MF

GOMF_tableSave <- enGO_HDR_merged_MF %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(p.adjust)) %>%
  dplyr::select(-plotpoint) %>%
  dplyr::mutate(geneID = gsub("/", ", ", geneID))

# Need to revert-back any simplifed Descritpions
# write.table(GOMF_tableSave, file = "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/final_plotsAndData/GOMF_HDRgenes.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Adding N2 genes, aliases, HOGs, and OGs
MF_N2_tableSave <- GOMF_tableSave %>%
  tidyr::separate_rows(geneID, sep = ",\\s") %>%
  dplyr::left_join(qn_table, by = c("geneID" = "QX1410")) %>%
  dplyr::left_join(n2_alias, by = c("N2" = "N2_WBGeneID")) %>%
  dplyr::arrange(Description, N2) %>%
  dplyr::left_join(n2_hogs_ogs, by = "N2") %>%
  dplyr::group_by(ID,Description,hdr_gene_term,hdr_gene_hit,bgd_term,bgd_hit,pvalue,p.adjust,qvalue,Count,EnrichRatio) %>%
  dplyr::summarise(
    geneID = paste(unique(geneID), collapse = ","),
    N2 = paste(unique(N2), collapse = ","),
    alias = paste(unique(alias), collapse = ","),
    HOGs = paste(unique(HOG), collapse = ","),
    OGs = paste(unique(OG), collapse = ","),
    .groups = "drop")

# write.table(MF_N2_tableSave, file = "/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/final_plotsAndData/GOMF_HDRgenes_updated_NOsulp7.tsv",  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


final_plot <- cowplot::plot_grid(
  plot_GO_BP, plot_GO_MF,
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("a", "b"),
  label_size = 14,
  label_fontface = "bold")
final_plot

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/processed_data/final_plotsAndData/GOBP_GOMF_Cbrig.png", final_plot,  width = 7.5, height = 7.5, dpi = 600)



########### No enrichment in non-HDR arm domains - need to combine with HDR arm genes and plot everything and have points different shapes for each respective set
# enGO_nHDR_merged <- clusterProfiler::enricher(
#   gene = nHD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "CC") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "CC") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX1410_merged_bgd_genes,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_nHDR_merged)
# 
# dotplot(enGO_nHDR_merged, showCategory = 40, title = "MF nHDRs")





































### Testing if all GO IDs are sourced from eggNOG-mapper or InterProScan
# Run eggNOG-mapper:
# emapper.py -m diamond -i QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.protein.fa \
# --data_dir /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/software/eggNOG_mapper/databases \
# --seed_ortholog_evalue 0.001 \
# --tax_scope eukaryota \
# --output_dir output \
# --temp_dir tmp \
# -o QX1410_longestIso \
# --cpu 48
# 
# grep -v "^#" QX1410_longestIso.emapper.annotations | awk -F'\t' -v OFS='\t' '$10 != "-" {print $1,$10}' > QX1410_longestIso.emapper.annotations.tsv

allQX <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/processed_data/QX1410_tranGene.tsv", col_names = c("tran","QX1410")) %>%
  dplyr::mutate(tran = paste0("transcript:",tran))

# egg <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/eggNOG/output/QX1410_longestIso.emapper.annotations.tsv", col_names = c("tran","GO")) %>%
#   tidyr::separate_rows(GO, sep = ',') %>%
#   dplyr::left_join(allQX, by = 'tran') %>%
#   dplyr::select(QX1410,GO) # 4,268

# QX_background <- unique(egg$QX1410) # XXX genes annotated with eggNOG-mapper
# 
# GO_annotations <- AnnotationDbi::select(GO.db, 
#                                          keys=unique(egg$GO),   
#                                          columns = c("TERM", "DEFINITION", "ONTOLOGY"),
#                                          keytype="GOID") %>% 
#   dplyr::rename(TERM = GOID, TERM_NAME = TERM) 
# 
# merged_ont <- egg %>%
#   dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) 
# 
# # BP
# enGO_HDR_merged <- clusterProfiler::enricher(
#   gene = HD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX_background,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_HDR_merged)
# 
# dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")
# 
# # MF
# enGO_HDR_merged_MF <- clusterProfiler::enricher(
#   gene = HD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX_background,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_HDR_merged_MF)
# 
# dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")


# 1e-4 e-value
# eggStringent <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/eggNOG/output/QX_1e-4.annotations.tsv", col_names = c("tran","GO")) %>%
#   tidyr::separate_rows(GO, sep = ',') %>%
#   dplyr::left_join(allQX, by = 'tran') %>%
#   dplyr::select(QX1410,GO)
# 
# QX_bckgrd_stringent <- unique(eggStringent$QX1410) # XXX genes annotated with eggNOG-mapper
# 
# GO_annotations <- AnnotationDbi::select(GO.db, 
#                                         keys=unique(eggStringent$GO),   
#                                         columns = c("TERM", "DEFINITION", "ONTOLOGY"),
#                                         keytype="GOID") %>% 
#   dplyr::rename(TERM = GOID, TERM_NAME = TERM) 
# 
# merged_ont <- eggStringent %>%
#   dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) 
# 
# # BP
# enGO_HDR_merged <- clusterProfiler::enricher(
#   gene = HD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX_bckgrd_stringent,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_HDR_merged)
# 
# dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")
# 
# # MF
# enGO_HDR_merged_MF <- clusterProfiler::enricher(
#   gene = HD_gene_vector,
#   TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
#   TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
#   universe = QX_bckgrd_stringent,
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05,
# )
# 
# head(enGO_HDR_merged_MF)
# 
# dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")




# 1e-4 e-value, but only arms as background
eggStringent <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/eggNOG/output/QX_1e-4.annotations.tsv", col_names = c("tran","GO")) %>%
  tidyr::separate_rows(GO, sep = ',') %>%
  dplyr::left_join(allQX, by = 'tran') %>%
  dplyr::select(QX1410,GO) %>% # 4,268
  dplyr::filter(QX1410 %in% arm_genes)

QX_bckgrd_stringent <- unique(eggStringent$QX1410) # only 1,255 genes...

GO_annotations <- AnnotationDbi::select(GO.db, 
                                        keys=unique(eggStringent$GO),   
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>% 
  dplyr::rename(TERM = GOID, TERM_NAME = TERM) 

merged_ont <- eggStringent %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) 

# BP
enGO_HDR_merged <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = QX_bckgrd_stringent,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged)

dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")

# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = QX_bckgrd_stringent,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")







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

ipr <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/briggsae/InterProScan/output/allApps_QX_IPR.tsv", col_names = c("tran","signature_accession","signature_description","IPR_accession","IPR_description","GO")) %>%
  dplyr::left_join(allQX, by = 'tran') %>%
  dplyr::select(-tran) %>%
  dplyr::select(QX1410,signature_accession,signature_description,IPR_accession,IPR_description,GO) 

# Unique InterPro per gene
ipr_gene <- ipr %>%
  dplyr::filter(!is.na(IPR_description) & IPR_description != "-") %>%
  dplyr::select(QX1410, IPR_accession, IPR_description) %>%
  dplyr::distinct(QX1410,IPR_accession, IPR_description) # 15,289 genes!

# Define universe & HDR membership (annotated-only universe) 
univ_genes <- unique(ipr_gene$QX1410)
hdr_genes  <- intersect(HD_gene_vector, univ_genes) # 3,427 genes 

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
  dplyr::arrange(FDR_p.adjust)

data_plt <- ipr_sig_gene_collapsed %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) +
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
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
        plot.margin = margin(t = 10, r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr

# data_plt <- ipr_sig_gene_collapsed %>%
#   dplyr::arrange(desc(enrich_ratio)) %>%
#   dplyr::slice_head(n = 20) %>%
#   dplyr::mutate(log10_p = -log10(FDR_p.adjust)) %>%
#   dplyr::arrange(enrich_ratio) %>%
#   dplyr::mutate(plotpoint = factor(row_number(), levels = rev(row_number()))) %>%
#   dplyr::arrange(desc(plotpoint)) %>%                                            
#   dplyr::mutate(IPR_description = factor(IPR_description, levels = IPR_description))
#                                    
# plot_ipr <- ggplot(data_plt, aes(x = enrich_ratio, y = IPR_description, fill = log10_p)) +
#   geom_col(width = 0.7) +
#   geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
#   scale_x_continuous(name = "Fold enrichment", expand = c(0.02, 0.02)) +
#   scale_fill_gradient(name = "p.adjust", labels = scales::label_number(accuracy = 20), low = "yellowgreen", high = "#53886C",) +
#   theme(
#     axis.text.y = element_text(size = 9, color = 'black'),
#     axis.text.x = element_text(size = 9, color = "black"),
#     axis.title.y = element_blank(),
#     legend.position = "right",
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_rect(color = 'black', fill = NA),
#     plot.title = element_text(size = 12, color = 'black')) +
#   labs(title = "Enriched IPR terms for genes in hyper-divergent regions")
# plot_ipr



##### NOW WITH JUST ARM GENES AS THE BACKGROUND
ipr_gene <- ipr %>%
  dplyr::filter(!is.na(IPR_description) & IPR_description != "-") %>%
  dplyr::select(QX1410, IPR_accession, IPR_description) %>%
  dplyr::distinct(QX1410,IPR_accession, IPR_description) %>%
  dplyr::filter(QX1410 %in% arm_genes) # 5,301 genes!

# Define universe & HDR membership (annotated-only universe) 
univ_genes <- unique(ipr_gene$QX1410)
hdr_genes  <- intersect(HD_gene_vector, univ_genes) # 3,427 genes 

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
  dplyr::arrange(FDR_p.adjust)

data_plt <- ipr_sig_gene_collapsed %>% dplyr::slice_head(n = 20) %>% dplyr::arrange(desc(FDR_p.adjust)) %>% dplyr::mutate(plotpoint = dplyr::row_number())

plot_ipr <- ggplot(data_plt) +
  geom_vline(xintercept = -log10(0.05), color='blue', linewidth=0.4) +
  geom_point(aes(x = -log10(FDR_p.adjust), y = plotpoint, size = enrich_ratio, fill = n_genes_HDR), shape = 21) +
  scale_y_continuous(breaks = data_plt$plotpoint, labels = data_plt$IPR_description, name = "", expand = c(0.02,0.02)) +
  scale_fill_gradient(low = "yellowgreen", high = "#53886C", breaks = c(round(min(data_plt$n_genes_HDR, na.rm = TRUE)), round((max(data_plt$n_genes_HDR, na.rm = TRUE) + min(data_plt$n_genes_HDR, na.rm = TRUE) ) / 2), round(max(data_plt$n_genes_HDR, na.rm = TRUE)))) +
  scale_size_continuous(range = c(0.5, 4), name = "Fold enrichment", breaks = pretty(data_plt$enrich_ratio, n = 3)) +
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
        plot.margin = margin(t = 10, r = 10, b = 10, l = 22, unit = "pt")) +
  # legend.spacing.x = unit(0.02, "cm"),
  # legend.key.size = unit(0.3, 'cm')) +
  guides(
    fill =  guide_colourbar(nrow=1, order = 1, title.position = "top", force = TRUE, barwidth = 5, barheight = 0.3),
    size = guide_legend(nrow=1, order = 2, title.position = "top", title.hjust = 1, force = TRUE)) +
  labs(title = "Enriched IPR terms for genes in HDRs",  x = expression(-log[10]~"(corrected p-value)"), size = "Fold enrichment", fill = "Gene count")
plot_ipr


# data_plt <- ipr_sig_gene_collapsed %>%
#   dplyr::arrange(desc(enrich_ratio)) %>%
#   dplyr::slice_head(n = 20) %>%
#   dplyr::mutate(log10_p = -log10(FDR_p.adjust)) %>%
#   dplyr::arrange(enrich_ratio) %>%
#   dplyr::mutate(plotpoint = factor(row_number(), levels = rev(row_number()))) %>%
#   dplyr::arrange(desc(plotpoint)) %>%                                            
#   dplyr::mutate(IPR_description = factor(IPR_description, levels = IPR_description))
# 
# plot_ipr <- ggplot(data_plt, aes(x = enrich_ratio, y = IPR_description, fill = log10_p)) +
#   geom_col(width = 0.7) +
#   geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
#   scale_x_continuous(name = "Fold enrichment", expand = c(0.02, 0.02)) +
#   scale_fill_gradient(name = "p.adjust", labels = scales::label_number(accuracy = 2), low = "yellowgreen", high = "#53886C",) +
#   theme(
#     axis.text.y = element_text(size = 9, color = 'black'),
#     axis.text.x = element_text(size = 9, color = "black"),
#     axis.title.y = element_blank(),
#     legend.position = "right",
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     panel.border = element_rect(color = 'black', fill = NA),
#     plot.title = element_text(size = 12, color = 'black')) +
#   labs(title = "Enriched IPR terms for genes in hyper-divergent regions")
# plot_ipr






##### GENE ONTOLOGY ENRICHMENT ANALYSIS FROM IPR GO IDs
go_ipr <- ipr %>%
  dplyr::filter(!is.na(GO) & GO != "-") %>%
  tidyr::separate_rows(GO, sep="\\|") %>%
  dplyr::filter(GO != "") %>%
  dplyr::distinct(QX1410, GO) %>% # 11,727 genes 
  dplyr::mutate(GO = str_remove_all(GO, "\\s*\\([^)]*\\)") |> str_squish())

IPR_GO_bckgrd <- unique(go_ipr$QX1410) # 11,727 genes

GO_annotations <- AnnotationDbi::select(GO.db, 
                                        keys=unique(go_ipr$GO),   
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>% 
  dplyr::rename(TERM = GOID, TERM_NAME = TERM) 

merged_ont <- go_ipr %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) 

# BP
enGO_HDR_merged <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged)

dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")

# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")


### Now with only arms as the background, not the entire genome
go_ipr_arms <- go_ipr %>% dplyr::filter(QX1410 %in% arm_genes)
IPR_GO_bckgrd_arms <- unique(go_ipr_arms$QX1410) # 3,900 genes

GO_annotations <- AnnotationDbi::select(GO.db, 
                                        keys=unique(go_ipr$GO),   
                                        columns = c("TERM", "DEFINITION", "ONTOLOGY"),
                                        keytype="GOID") %>% 
  dplyr::rename(TERM = GOID, TERM_NAME = TERM) 

merged_ont <- go_ipr %>%
  dplyr::left_join(GO_annotations, by = c("GO" = "TERM")) 

# BP
enGO_HDR_merged <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "BP") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged)

dotplot(enGO_HDR_merged, showCategory = 40, title = "BP HDRs")

# MF
enGO_HDR_merged_MF <- clusterProfiler::enricher(
  gene = HD_gene_vector,
  TERM2GENE = merged_ont %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(GO,QX1410),
  TERM2NAME = GO_terms_merged %>% dplyr::filter(ONTOLOGY == "MF") %>% dplyr::select(TERM,TERM_NAME),
  universe = IPR_GO_bckgrd_arms,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
)

head(enGO_HDR_merged_MF)

dotplot(enGO_HDR_merged_MF, showCategory = 40, title = "MF HDRs")

