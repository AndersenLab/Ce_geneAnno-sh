library(dplyr)
library(tidyr)
library(stringr)
library(ape)

gff <- ape::read.gff("/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/csq/c_elegans.PRJNA13758.WS283.csq.gff3")


PC_genes <- gff %>%
  dplyr::filter(grepl("protein_coding", attributes) & type == "gene") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  dplyr::mutate(L1_ID=gsub("ID=","",pre)) %>%
  dplyr::select(-pre,-post) 

PC_mRNA <- gff %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  dplyr::select(-post) %>%
  tidyr::separate(pre,into=c("L2_ID","L1_ID"),sep=";") %>%
  dplyr::mutate(L2_ID=gsub("ID=","",L2_ID)) %>%
  dplyr::mutate(L1_ID=gsub("Parent=","",L1_ID)) %>%
  dplyr::filter(L1_ID %in% PC_genes$L1_ID)

L3_ofeat <- c("exon","three_prime_UTR","five_prime_UTR")

PC_L3 <- gff %>%
  dplyr::filter(type %in% L3_ofeat) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","",attributes)) %>%
  dplyr::filter(L2_ID %in% PC_mRNA$L2_ID)

PC_CDS <- gff %>%
  dplyr::filter(type=="CDS") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  tidyr::separate(pre,into=c("L3_ID","L2_ID"),sep=";") %>%
  dplyr::select(-post,-L3_ID) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","",L2_ID)) %>%
  dplyr::filter(L2_ID %in% PC_mRNA$L2_ID)

PC_gff <- rbind(
  PC_genes %>% dplyr::select(-L1_ID),
  PC_mRNA %>% dplyr::select(-L1_ID,-L2_ID),
  PC_L3 %>% dplyr::select(-L2_ID),
  PC_CDS %>% dplyr::select(-L2_ID)) %>% 
  dplyr::arrange(seqid,start) %>%
  dplyr::mutate(score=".") %>%
  dplyr::mutate(phase=as.character(phase)) %>%
  dplyr::mutate(phase=ifelse(is.na(phase),".",phase))


write.table(PC_gff,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.gff3", quote = F, row.names = F, col.names = F, sep = '\t')
