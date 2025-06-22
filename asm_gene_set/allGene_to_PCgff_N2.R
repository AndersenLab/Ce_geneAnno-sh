library(dplyr)
library(tidyr)
library(stringr)
library(ape)

# N2
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


# QX1410 
gffb <- ape::read.gff("/vast/eande106/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/csq/QX1410.update.April2025.noWBGeneID.csq.gff3")

PC_genesb <- gffb %>%
  dplyr::filter(grepl("protein_coding", attributes) & type == "gene") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  dplyr::select(-pre,-post) %>%
  tidyr::separate(attributes,into=c("L1_ID","post"),sep=";",remove = F) %>%
  dplyr::select(-post) %>%
  dplyr::mutate(L1_ID = gsub("ID=","", L1_ID)) 

PC_mRNAb <- gffb %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep="ID=",remove = F) %>%
  dplyr::select(-pre) %>%
  tidyr::separate(post,into=c("L2_ID","L1_ID"),sep=";") %>%
  dplyr::mutate(L1_ID=gsub("Parent=","",L1_ID)) %>%
  dplyr::filter(L1_ID %in% PC_genesb$L1_ID)

L3_ofeatb <- c("exon","three_prime_UTR","five_prime_UTR")

PC_L3b <- gffb %>%
  dplyr::filter(type %in% L3_ofeat) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","", attributes)) %>%
  tidyr::separate(L2_ID,into=c("L2_ID","post"),sep=";",remove = F) %>%
  dplyr::select(-L2_ID) %>%
  dplyr::filter(post %in% PC_mRNAb$L2_ID)

PC_CDSb <- gffb %>%
  dplyr::filter(type=="CDS") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  tidyr::separate(pre,into=c("L3_ID","L2_ID"),sep=";") %>%
  dplyr::select(-post,-L3_ID) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","",L2_ID)) %>%
  dplyr::filter(L2_ID %in% PC_mRNAb$L2_ID)

PC_gffb <- rbind(
  PC_genesb %>% dplyr::select(-L1_ID),
  PC_mRNAb %>% dplyr::select(-L1_ID,-L2_ID),
  PC_L3b %>% dplyr::select(-post),
  PC_CDSb %>% dplyr::select(-L2_ID)) %>% 
  dplyr::arrange(seqid,start) %>%
  dplyr::mutate(score=".") %>%
  dplyr::mutate(phase=as.character(phase)) %>%
  dplyr::mutate(phase=ifelse(is.na(phase),".",phase))

write.table(PC_gffb,"/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/QX1410.update.April2025.noWBGeneID.csq.ONLYPC.gff3", quote = F, row.names = F, col.names = F, sep = '\t')


# AF16
gffa <- ape::read.gff("/vast/eande106/data/c_briggsae/genomes/PRJNA10731/WS276/csq/c_briggsae.PRJNA10731.WS276.csq.gff3.gz")

PC_genesa <- gffa %>%
  dplyr::filter(grepl("protein_coding", attributes) & type == "gene") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  dplyr::mutate(L1_ID=gsub("ID=","",pre)) %>%
  dplyr::select(-pre,-post) 

PC_mRNAa <- gffa %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  dplyr::select(-post) %>%
  tidyr::separate(pre,into=c("L2_ID","L1_ID"),sep=";") %>%
  dplyr::mutate(L2_ID=gsub("ID=","",L2_ID)) %>%
  dplyr::mutate(L1_ID=gsub("Parent=","",L1_ID)) %>%
  dplyr::filter(L1_ID %in% PC_genesa$L1_ID) %>%
  dplyr::mutate(L2_ID = gsub("CDS:","",L2_ID))

L3_ofeat <- c("exon","three_prime_UTR","five_prime_UTR")

PC_L3a <- gffa %>%
  dplyr::filter(type %in% L3_ofeat) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","",attributes)) %>%
  dplyr::mutate(L2_ID = gsub("CDS:","",L2_ID)) %>%
  dplyr::filter(L2_ID %in% PC_mRNAa$L2_ID)

PC_CDSa <- gffa %>%
  dplyr::filter(type=="CDS") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F)%>%
  tidyr::separate(pre,into=c("L3_ID","L2_ID"),sep=";") %>%
  dplyr::select(-post,-L3_ID) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","",L2_ID)) %>%
  dplyr::filter(L2_ID %in% PC_mRNAa$L2_ID)

PC_gffa <- rbind(
  PC_genesa %>% dplyr::select(-L1_ID),
  PC_mRNAa %>% dplyr::select(-L1_ID,-L2_ID),
  PC_L3a %>% dplyr::select(-L2_ID),
  PC_CDSa %>% dplyr::select(-L2_ID)) %>% 
  dplyr::arrange(seqid,start) %>%
  dplyr::mutate(score=".") %>%
  dplyr::mutate(phase=as.character(phase)) %>%
  dplyr::mutate(phase=ifelse(is.na(phase),".",phase))


write.table(PC_gffa,"/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/c_briggsae.PRJNA10731.WS276.csq.ONLYPC.gff3", quote = F, row.names = F, col.names = F, sep = '\t')




# NIC58 
gfft <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/tropicalis/prep/NIC58.update.April2025.noWBGeneID.csq.gff3")

PC_genest <- gfft %>%
  dplyr::filter(grepl("protein_coding", attributes) & type == "gene") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  dplyr::select(-pre,-post) %>%
  tidyr::separate(attributes,into=c("L1_ID","post"),sep=";",remove = F) %>%
  dplyr::select(-post) %>%
  dplyr::mutate(L1_ID = gsub("ID=","", L1_ID)) 

PC_mRNAt <- gfft %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep="ID=",remove = F) %>%
  dplyr::select(-pre) %>%
  tidyr::separate(post,into=c("L2_ID","L1_ID"),sep=";") %>%
  dplyr::mutate(L1_ID=gsub("Parent=","",L1_ID)) %>%
  dplyr::filter(L1_ID %in% PC_genest$L1_ID)

L3_ofeatt <- c("exon","three_prime_UTR","five_prime_UTR")

PC_L3t <- gfft %>%
  dplyr::filter(type %in% L3_ofeatt) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","", attributes)) %>%
  tidyr::separate(L2_ID,into=c("L2_ID","post"),sep=";",remove = F) %>%
  dplyr::select(-L2_ID) %>%
  dplyr::filter(post %in% PC_mRNAt$L2_ID)

PC_CDSt <- gfft %>%
  dplyr::filter(type=="CDS") %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=";Name=",remove = F) %>%
  tidyr::separate(pre,into=c("L3_ID","L2_ID"),sep=";") %>%
  dplyr::select(-post,-L3_ID) %>%
  dplyr::mutate(L2_ID=gsub("Parent=","",L2_ID)) %>%
  dplyr::filter(L2_ID %in% PC_mRNAt$L2_ID)

PC_gfft <- rbind(
  PC_genest %>% dplyr::select(-L1_ID),
  PC_mRNAt %>% dplyr::select(-L1_ID,-L2_ID),
  PC_L3t %>% dplyr::select(-post),
  PC_CDSt %>% dplyr::select(-L2_ID)) %>% 
  dplyr::arrange(seqid,start) %>%
  dplyr::mutate(score=".") %>%
  dplyr::mutate(phase=as.character(phase)) %>%
  dplyr::mutate(phase=ifelse(is.na(phase),".",phase))

write.table(PC_gfft,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/GO_enrichment/tropicalis/prep/NIC58.update.April2025.noWBGeneID.csq.ONLYPC.gff3", quote = F, row.names = F, col.names = F, sep = '\t')







### LOAD IN LONGEST-ISOFORM OUTPUT FROM AGAT ###
# ^^^^^^^^^^^^^^ NOT NEEDED!!!! AGAT IS FINE - it keeps multiple transcripts if the mRNA is the same length, but one entry won't have CDSs but the other does... so when you 
# translate it only translates one of the transcripts
# longestIso_gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3")
# 
# 
# grouping <- longestIso_gff %>%
#   dplyr::mutate(id = str_extract(attributes, "WBGene[0-9]+"), tran_id = str_extract(attributes, "transcript:[^;]+") %>% str_remove("transcript:")) %>%
#   dplyr::mutate(mRNA_len = ifelse(type == "mRNA", (end - start), NA))
# 
# longest_transcripts <- grouping %>%
#   dplyr::filter(type == "mRNA") %>%
#   dplyr::group_by(id) %>%
#   dplyr::slice_max(order_by = mRNA_len, n = 1, with_ties = FALSE) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(id, tran_id)
# 
# longest_gff <- grouping %>%
#   dplyr::filter(tran_id %in% longest_transcripts$tran_id | (type == "gene" & id %in% longest_transcripts$id)) %>%
#   dplyr::mutate(score=".") %>%
#   dplyr::mutate(phase=as.character(phase)) %>%
#   dplyr::mutate(phase=ifelse(is.na(phase),".", phase)) %>%
#   dplyr::select(-id, -tran_id, -mRNA_len)
# 
# 
# write.table(longest_gff,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3", quote = F, row.names = F, col.names = F, sep = '\t')
