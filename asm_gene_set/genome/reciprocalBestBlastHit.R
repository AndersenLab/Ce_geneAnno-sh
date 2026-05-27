library(stringr)
library(dplyr)
library(tidyr)
library(ape)
library(readr)

# Read in GFFs and extract transcript name and gene name pairs
QXgff <- ape::read.gff("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/prep/QX1410.update.April2025.noWBGeneID.csq.ONLYPC.longest.gff3")
QXgff <- QXgff %>%
  dplyr::filter(type=="mRNA") %>%
  tidyr::separate(attributes,into=c("pre","post"), sep="ID=",remove = F) %>%
  dplyr::select(-pre) %>%
  tidyr::separate(post,into=c("L2_ID","L1_ID"),sep=";") %>%
  dplyr::mutate(L1_ID=gsub("Parent=","",L1_ID)) %>%
  dplyr::rename(QXtran = L2_ID, QXgene = L1_ID) %>%
  dplyr::select(QXtran, QXgene)

N2gff <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/longest_isoform/c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3")
N2gff <- N2gff %>% 
  dplyr::filter(type=="mRNA") %>% 
  tidyr::separate(attributes,into=c("N2tran","N2gene","Name","Other"),sep=";",extra="merge") %>%
  dplyr::mutate(N2tran = gsub(".*=","", N2tran)) %>%
  dplyr::mutate(N2gene =gsub(".*=","", N2gene)) %>%
  dplyr::select(N2tran, N2gene)


# Read in QX1410 against N2 BLASTP output
QX_blastp <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/RBBH/QXagainstN2.recipro.pb.out", 
                            col_names = c("QXtran", "N2tran", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "eval", "bitscore"))

QX_bestHit <- QX_blastp %>%
  dplyr::filter(eval < 1e-5) %>%
  dplyr::left_join(QXgff, by = "QXtran") %>%
  dplyr::group_by(QXtran) %>%
  dplyr::arrange(eval,  desc(bitscore), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(N2gff, by = "N2tran") %>%
  dplyr::select(N2gene, QXgene)


# Read in N2 against QX1410 BLASTP output
N2_blastp <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/RBBH/N2againstQX.pb.out",
                             col_names = c("N2tran", "QXtran", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "eval", "bitscore"))

N2_bestHit <- N2_blastp %>%
  dplyr::filter(eval < 1e-5) %>%
  dplyr::left_join(N2gff, by = "N2tran") %>%
  dplyr::group_by(N2gene) %>%
  dplyr::arrange(eval, desc(bitscore), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(QXgff, by = "QXtran") %>%
  dplyr::select(N2gene, QXgene)


RBBH <- dplyr::inner_join(QX_bestHit, N2_bestHit, by = c("N2gene", "QXgene")) %>% 
  dplyr::mutate(across(everything(), ~ gsub("gene:", "", .)))

write.table(RBBH, "/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/briggsae/RBBH/N2_QX1410_RBBH.tsv", sep="\t", col.names=FALSE, row.names = FALSE, quote=F)
