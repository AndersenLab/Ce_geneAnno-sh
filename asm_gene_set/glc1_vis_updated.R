library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(grid)



##### This section of code is taken from Nic's old haplotypePlotter R script ####
setwd("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/glc1_variation/HDR_haplotypePlotter")

#read collapsed reference HDRs
collapsed_ff <- readr::read_tsv("./input/HDR_5kbclust_collapsed_wFreq.tsv")

#read strain-specific HDRs
# all_SR_calls <- readr::read_tsv("./input/HDR_allStrain_5kbclust_1IBfilt.tsv")

#read all pairwise genome coordinate comparisons
### This will need to eventually change to all pairwise alignments among all WSs to the reference
transformed_coords <- readr::read_tsv("/vast/eande106/projects/Nicolas/c.elegans/reference_genealn/N2vCB/N2_hifi_transformed2.tsv",col_names = F) %>% dplyr::mutate(STRAIN="CB4856")

colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")

#read concatentated gene models of every genome
gffCat <- ape::read.gff("./input/all_LRiso.gff")
gffCat1 <- ape::read.gff("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/braker_runs/gff/CB4856.braker.gff3") %>% dplyr::mutate(STRAIN="CB4856")
gffCat2 <- ape::read.gff("/vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3") %>% dplyr::mutate(STRAIN="N2")
gffCat <- rbind(gffCat1,gffCat2)

#read ortholog relationships among gene models
#orthos <- readr::read_tsv("./input/Orthogroups.tsv")
orthos <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/N2xCB/OrthoFinder/Results_Mar17/Orthogroups/Orthogroups.tsv")
strainCol <- colnames(orthos)
strainCol_c1 <- gsub(".braker.protein","",strainCol)
strainCol_c2 <- gsub("_WS283.protein","",strainCol_c1)


#set your target  coordinates
#GLC-1
hdr_chrom = "V"
hdr_start_pos = 16115967
hdr_end_pos = 16276907

#offset lets you explore adjacent regions
offset = 0
hap_chrom = hdr_chrom
hap_start = hdr_start_pos - offset
hap_end = hdr_end_pos + offset #+ 10000

#use reference coordinates from g2g alginments to pull the contigs that contain the alt haplotypes for the HDR
hap_coords <- transformed_coords %>%
  dplyr::filter((REF == hap_chrom & hap_start >= S1 & hap_start <= E1 ) | 
                  (REF == hap_chrom & hap_end >= S1 & hap_end <= E1) | 
                  (REF == hap_chrom & S1 >= hap_start & E1 <= hap_end)) %>%
  dplyr::mutate(inv=ifelse(S2>E2,T,F)) 
# dplyr::select(-S2,-E2) %>%
# dplyr::rename(S2=newS2,E2=newE2)

# naive visualization of g2g alignments for the target region
# multiple contigs may map to the REF region, we need to filter those!
ggplot(hap_coords) + geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=inv)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)")

#keep only the contig with the largest extent of alignment with the REF HDR
tigFilt <- hap_coords %>% 
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(HIFI) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L2)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::ungroup()

#important diagnostic plot
#after the optimal contig is selected, we can see if the HDR boundaries are within the alignment
#for proper visualization, the HDR (grey box in plot) needs to be encompassed by the selected contig
#otherwise some haplotypes will be truncated 
#a step to drop genomes with incomplete coverage of the HDR could be added
#SEA-1 locus behaves well
#have in mind that the alignments look fragmented because of sequence divergence, but the contig is linear in the genome file
ggplot(tigFilt) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA))

#get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
HV_boundary <- tigFilt %>%
  dplyr::mutate(refStart=min(S1),refEnd=max(E1)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(boundStart=min(S2), boundEnd=max(E2)) %>%
  dplyr::distinct(STRAIN, .keep_all = T) %>%
  dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd) %>%
  dplyr::ungroup()

#filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
filtGff <- gffCat %>%
  dplyr::filter(type=="gene") %>%
  dplyr::filter(seqid %in% HV_boundary$HIFI) %>%
  dplyr::left_join(HV_boundary,by=c('seqid'='HIFI',"STRAIN")) %>%
  dplyr::filter((start >= boundStart & start <= boundEnd) | 
                  (end >= boundStart & end <= boundEnd)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ngene=n()) %>%
  dplyr::arrange(ngene) %>%
  dplyr::mutate(gid=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(boundSize=abs(boundStart-boundEnd))

#extract the REF genes 
N2Start = filtGff$refStart[1]
N2End = filtGff$refEnd[1]
N2_genes <- gffCat %>%
  dplyr::filter(seqid==hap_chrom & type=="gene") %>%
  dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | 
                  (end >= hap_start & end <= hap_end)) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  dplyr::mutate(seqname=paste0("Transcript_",seqname)) %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name))  

#extract the REF protein-coding transcripts
N2_tran <- gffCat %>%
  dplyr::filter(seqid==hap_chrom & type=="mRNA") %>%
  tidyr::separate(attributes, into=c("ID","Parent","Name","wormpep","locus"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=Transcript:","",ID)) %>%
  tidyr::separate(ID,into = c("fosmid","tseqID",'tranum'),sep="\\.",remove = F) %>%
  dplyr::mutate(tseqname=paste0("Transcript_",fosmid,".",tseqID,".",tranum)) %>%
  dplyr::mutate(Parent=gsub("Parent=Gene:","",Parent)) %>%
  dplyr::filter(Parent %in% N2_genes$Name) %>%
  dplyr::select(tseqname,Parent) %>%
  dplyr::rename(tranname=tseqname) %>%
  dplyr::left_join(N2_genes,by=c('Parent'='Name'))

# #get gene list
HV_genelist <- N2_genes$seqname

#get alt gene names/aliases
aliases <- N2_genes %>%
  dplyr::select(seqname,post2) %>%
  tidyr::separate(post2,into=c("rem","aliases"),sep=";Alias=") %>%
  dplyr::select(-rem) %>%
  tidyr::separate(aliases,into=c('locus_name',"alias2"),sep=',') %>%
  dplyr::select(-alias2)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
ggplot(N2_genes) + geom_rect(aes(xmin=start,xmax=end,ymin=1,ymax=2))

#filter orthologous groups using REF genes
#this will establish your orthology relationships between REF and WILD haplotypes
filtOrthos <- orthos %>%
  dplyr::filter(grepl(paste(HV_genelist,collapse="|"),N2_WS283.protein)) %>%
  dplyr::mutate(N2 = strsplit(as.character(N2_WS283.protein), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end),by=c("N2"="tranname")) %>%
  dplyr::filter(!is.na(seqid)) %>%
  dplyr::distinct(Orthogroup,.keep_all = T) %>%
  dplyr::mutate(N2=seqname)

#generate a lookup table (all_ortho_pairs) which contains all pairwise gene orthologs between REF and WILD
orthoList <- list()
for (i in 2:length(strainCol)) {
  tmp <- filtOrthos %>%
    dplyr::select(N2,strainCol[i]) %>%
    dplyr::rename(tmpSel=strainCol[i]) %>%
    dplyr::mutate(newSel = strsplit(as.character(tmpSel), ",")) %>%
    tidyr::unnest(newSel) %>%
    dplyr::mutate(newSel=trimws(newSel)) %>%
    dplyr::select(-tmpSel) %>%
    dplyr::mutate(newSel=paste0(strainCol_c2[i],".",newSel)) %>%
    tidyr::separate(newSel,into = c("STRAIN","geneid",'tranid'),sep="\\.") %>%
    dplyr::mutate(N2=gsub("([a-zA-Z])$", "",N2)) ###### ADD TO HAP PLOTTER ##########
  
  orthoList[[i-1]] <- tmp
}
all_ortho_pairs  <- ldply(orthoList,data.frame) %>%
  dplyr::mutate(geneid=ifelse(STRAIN=='N2',paste0(geneid,".",tranid),geneid))

#prep GFF fields for join with ortholog table
gffFilt_ortho <- gffCat %>%
  dplyr::filter((!source=="WormBase") & type=="gene") %>%
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::mutate(attributes=gsub(";","",attributes)) %>%
  dplyr::rename(geneid=attributes) %>%
  dplyr::mutate(tigID1=seqid)
#%>%
  #tidyr::separate(seqid,into = c("STRAIN","tigID1",'tigID2'),sep = '\\.',remove = F)

# #join ortholog table to GFF coordinates
#seqid delineates the contig (thus the wild strain), which sets the Y position in the haplotype plot
#keep the longest transcript per gene
orthoPairs_wcoord <- all_ortho_pairs %>%
  dplyr::left_join(gffFilt_ortho,by = c("STRAIN","geneid")) %>%
  dplyr::arrange(seqid) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(ypos=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(gsize==max(gsize)) %>%
  dplyr::ungroup()

#similar to DF above, but
#some contigs are inverted, and gene models need to be adjusted for orientation
#some contigs have partial alignments that contain genes further away from the local HDR
#using lead() we can remove genes that have large (>5e4) jumps in position
#this trimming process may be the cause of misbehavior in plots of other HDRs
orthoPairs_wcoord_trimmer <- all_ortho_pairs %>%
  dplyr::left_join(gffFilt_ortho,by = c("STRAIN","geneid")) %>%
  dplyr::arrange(seqid) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(ypos=cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(gsize==max(gsize)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(hap_coords %>%
                     dplyr::select(HIFI,inv) %>%
                     dplyr::distinct(HIFI,.keep_all=T),by=c("seqid"="HIFI")) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(newStart=ifelse(inv==T,abs(end-(max(end))),start),newEnd=ifelse(inv==T,abs(start-(max(end))),end)) %>%
  #dplyr::select(-start,-end) %>%
  #dplyr::rename(start=newStart,end=newEnd) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(seqid,newStart) %>%
  dplyr::filter(!(is.na(STRAIN))) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDist=abs(newEnd-lead(newStart))) %>%
  dplyr::mutate(trm=ifelse(leadDist>5e4,"S",NA)) %>%
  #dplyr::mutate(allPass=ifelse(any(trm=="S"),"FIX","KEEP")) %>%
  dplyr::mutate(check=ifelse(!(is.na(trm)),row_number(),NA)) %>%
  dplyr::mutate(fill=ifelse(STRAIN=="N2",0, check)) %>%
  tidyr::fill(fill,.direction = 'up' ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,fill) %>%
  dplyr::mutate(gsize_filt=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(filt=ifelse(gsize_filt==max(gsize_filt),F,T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(filt==F)


#dplyr::mutate(fill=ifelse(STRAIN=="N2",0,fill_run(check, run_for_first = T)))#this is the problematic line
# dplyr::filter(is.na(trm) | trm=="NR")

#pull the boundaries of the ortholgous genes
pullBound <- orthoPairs_wcoord_trimmer %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(boundStart = min(start), boundEnd = max(end)) %>%
  dplyr::distinct(seqid,.keep_all = T) %>%
  dplyr::select(STRAIN,seqid,boundStart,boundEnd) %>%
  dplyr::ungroup()

#find the bound genes for each strain that are not orthologous
boundGenes <- gffCat %>%
  dplyr::filter(type=="gene" & seqid %in% pullBound$seqid) %>%
  dplyr::left_join(pullBound,by=c("seqid","STRAIN")) %>%
  dplyr::group_by(seqid) %>%
  dplyr::filter(start >= boundStart & start <= boundEnd) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(geneid=gsub("ID=","",attributes)) %>%
  dplyr::mutate(geneid=gsub(";","",geneid)) %>%
  dplyr::left_join(orthoPairs_wcoord,by=c('geneid','STRAIN')) %>%
  dplyr::filter(is.na(N2)) %>%
  dplyr::select(seqid.x, start.x, end.x, geneid, strand.x) %>%
  dplyr::rename(seqid=seqid.x,start=start.x,end=end.x,strand=strand.x) %>%
  dplyr::mutate(N2="non-ortho")


#bind ortholgous and non orthologous genes for each strain, and transform the coordinates of genes to a midle axis center
plotCoords <- rbind(boundGenes %>% dplyr::mutate(geneid=gsub(";","",geneid)),orthoPairs_wcoord_trimmer %>% dplyr::select(seqid,start,end,geneid,strand,N2)) %>%
  dplyr::mutate(ortho_status=ifelse(N2=='non-ortho',F,T)) %>%
  dplyr::left_join(hap_coords %>%
                     dplyr::select(HIFI,inv) %>%
                     dplyr::distinct(HIFI,.keep_all=T),by=c("seqid"="HIFI")) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(newStart=ifelse(inv==T,abs(end-(max(end))),start),newEnd=ifelse(inv==T,abs(start-(max(end))),end)) %>%
  dplyr::select(-start,-end) %>%
  dplyr::rename(start=newStart,end=newEnd) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(seqid,start) %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(leadDist=abs(end-lead(start))) %>%
  dplyr::mutate(trm=ifelse(leadDist>5e4 | lag(leadDist) >5e4,"R","NR")) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.na(trm) | trm=="NR") %>%
  dplyr::group_by(seqid) %>%
  dplyr::mutate(newstart=start-min(start),newend=end-min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(STRAIN="CB4856") %>%
  dplyr::mutate(tig1=seqid) %>%
  #tidyr::separate(seqid,into = c("STRAIN","tig1","tig2"),sep = "\\.",remove = F) %>%
  dplyr::select(-tig1,-inv,-start,-end,-trm,-leadDist) %>%
  dplyr::select(seqid,STRAIN,N2,ortho_status,newstart,newend,strand)

#get N2 genes
N2ad <- N2_genes %>%
  dplyr::mutate(STRAIN="N2",ortho_status=T) %>%
  dplyr::rename(N2=seqname) %>%
  dplyr::mutate(shift=min(start)) %>%
  dplyr::mutate(newstart=start-shift,newend=end-shift) %>%
  dplyr::select(seqid,STRAIN,N2,ortho_status,newstart,newend,shift,strand) %>%
  dplyr::mutate(sp = ifelse(strand == "+", 2.25,1.75)) 











######## Here is where the TSVs you create in bash will be directly loaded in and used ######## 
# Load in SV and SNV calls by paftools
df = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/elegans.merged.1kbCOV.5kbALIGN.annotated.final.vcf")

CB4856var <- df %>%
  dplyr::rename(CHROM = '#CHROM', start = POS) %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856) %>%
  dplyr::filter(CHROM == "V", (start >= 16115967 & start <= 16276907)) %>%
  dplyr::filter(CB4856 != "./.") 

plot_df <- CB4856var %>%
  dplyr::mutate(
    lenDEL = case_when(INFO == "DEL" ~ (end - start), TRUE ~ NA_real_),
    lenINS = case_when(INFO == "INS" ~ (end - start), TRUE ~ NA_real_))

# SNPs <- plot_df %>%
#   dplyr::filter(INFO == "SNP") %>%
#   dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856) %>%
#   dplyr::filter(!grepl(",", ALT))

deletions <- plot_df %>%
  dplyr::filter(INFO == "DEL") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenDEL) %>%
  dplyr::filter(lenDEL >= 50)

insertions <- plot_df %>%
  dplyr::filter(INFO == "INS") %>%
  dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenINS) %>%
  dplyr::filter(lenINS >= 50)



# Load in SNPs called by GATK pipeline with SR data
SR = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/SR.GATK.CB4856.final.vcf")

GATK <- SR %>%
  dplyr::filter(POS >= 16115967 & POS <= 16276907) %>%
  dplyr::rename(start = POS) %>%
  dplyr::mutate(end = start+1) %>%
  dplyr::filter(CB4856 != './.' & CB4856 != '0/0')

indels = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/SR.GATK.CB4856.indels.final.vcf")

G_del <- indels %>%
  dplyr::filter(START >= 16115967 & START <= 16276907) %>%
  dplyr::filter(ANNO == "DEL") %>%
  dplyr::filter() %>%
  dplyr::filter((END - START) < 50)

G_ins <- indels %>%
  dplyr::filter(START >= 16115967 & START <= 16276907) %>%
  dplyr::filter(ANNO == "INS") %>%
  dplyr::filter((END - START) < 50)

hist_data <- ggplot_build(
  ggplot(GATK, aes(x = start / 1e6)) + geom_histogram(binwidth = 0.0002)
)$data[[1]]

y_max <- max(hist_data$count)
print(y_max)
y_scaler <- 0.2 / y_max

hist <- ggplot(GATK, aes(x = start / 1e6)) +
  geom_histogram(binwidth = 0.0002, aes(y = ..count.. * y_scaler, fill = 'SNVs (GATK)'), alpha = 0.8) + 
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  scale_fill_manual(
    name = "",
    values = c("SNVs (GATK)" = "purple"),
    breaks = c("SNPs (GATK)")
  ) +
  scale_y_continuous(limits = c(0, 0.2)) 

hist


plot1 <- ggplot() +
  geom_segment(data = N2ad, aes(x = (hdr_start_pos)/1e6, xend = (hdr_end_pos)/1e6, y = 2.25, yend = 2.25), 
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")) + 
  geom_segment(data = N2ad, aes(x = (hdr_end_pos)/1e6, xend = (hdr_start_pos)/1e6 - 0.002, y = 1.75, yend = 1.75), 
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")) + 
  geom_rect(
    data = N2ad,
    aes(xmin = (newstart + shift)/1e6, xmax = (newend + shift)/1e6, ymin = sp + 0.2, ymax = sp - 0.2, fill = N2),
    color = 'black'
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.background = element_blank()
  )

plot1

plot2 <- ggplot() +
  # geom_rect(data = plot_df, aes(xmin = min(start)/1e6, xmax = max(start)/1e6, ymin = 1.75, ymax = 2.25, fill = 'QTL')) +
  # geom_rect(data = plot_df, aes(xmin = 16219634/1e6, xmax = 16221917/1e6, ymin = 1.75, ymax = 2.25, fill = 'glc-1')) +
  
  geom_rect(data = hist_data, aes(xmin = xmin - 0.0001, xmax = xmin + 0.0001, ymin = 1.5,  ymax = count * y_scaler + 1.5, fill = 'SNVs (GATK)'), alpha = 0.8) +
  
  geom_rect(data = GATK, aes(xmin = start/1e6 - 0.00001, xmax = end/1e6, ymin = 1.0, ymax = 1.5, fill = 'SNVs (GATK)')) + 
  geom_rect(data = deletions, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.25, ymax = 0.75, fill = 'Deletions')) +
  geom_rect(data = G_del, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.25, ymax = 0.75, fill = 'Deletions')) +
  geom_rect(data = G_ins, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.25, ymax = 0.75, fill = 'Insertions')) +
  geom_rect(data = insertions, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.25, ymax = 0.75, fill = 'Insertions')) +
  
  scale_x_continuous(name = "Genomic Position (Mb)", labels = scales::number_format(scale = 1, accuracy = 0.01)) +

  scale_fill_manual(
    name = "",
    values = c("QTL" = "darkolivegreen", 'glc-1' = 'darkolivegreen3', "SNVs (GATK)" = "purple", "Insertions" = "blue", "Deletions" = "red"),
    breaks = c("QTL", "glc-1", "SNPs (GATK)", "Insertions", "Deletions")
  ) + 
  
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line.x = element_line(),
    legend.position = "none",
    panel.background = element_blank()
  )
  
plot2


combined_plot <- plot_grid(plot1, plot2, ncol = 1, align = "v", rel_heights = c(0.4, 2))

combined_plot


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc_1_wGenes.jpg", combined_plot, dpi=600, width = 7.5, height = 6)



N2ad$sp <- as.factor(N2ad$sp)

plot3 <- ggplot() +
  geom_segment(data = N2ad, aes(x = (hdr_start_pos)/1e6, xend = (hdr_end_pos)/1e6, y = 2.25, yend = 2.25)) +
  geom_rect(
    data = N2ad,
    aes(xmin = (newstart + shift)/1e6, xmax = (newend + shift)/1e6, ymin = 2.225, ymax = 2.275, fill = sp),
    color = 'black'
  ) +
  scale_fill_manual(values = c("gray50", "gray90")) + # Watson is gray90 (lightgray)
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none",
    panel.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

plot3

plot4 <- ggplot() +
  geom_rect(data = GATK, aes(xmin = start/1e6 - 0.00001, xmax = end/1e6, ymin = 0.6, ymax = 0.85, fill = 'SNVs (GATK)')) + 
  geom_rect(data = deletions, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.25, ymax = 0.5, fill = 'Deletions')) +
  geom_rect(data = G_del, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.25, ymax = 0.5, fill = 'Deletions')) +
  geom_rect(data = G_ins, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.25, ymax = 0.5, fill = 'Insertions')) +
  geom_rect(data = insertions, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.25, ymax = 0.5, fill = 'Insertions')) +
  
  scale_x_continuous(name = "Genomic position (Mb)", labels = scales::number_format(scale = 1, accuracy = 0.01)) +
  
  scale_fill_manual(
    name = "",
    values = c("QTL" = "darkolivegreen", 'glc-1' = 'darkolivegreen3', "SNVs (GATK)" = "black", "Insertions" = "blue", "Deletions" = "red"),
    breaks = c("QTL", "glc-1", "SNPs (GATK)", "Insertions", "Deletions")
  ) + 
  
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line.x = element_line(),
    legend.position = "none",
    panel.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

plot4


combined_plot2 <- plot_grid(plot3, plot4, ncol = 1, align = "v", rel_heights = c(0.4, 0.6))

combined_plot2

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc_1_wGenes_updated.jpg", combined_plot2, dpi=600, width = 7.5, height = 1.8)

# gfirst <- "N2.F44G3.11"
# glast <- "T13F3.7"
# 
# #store the N2 coordinate shift
# N2_shift <- c("N2",unique(N2ad$shift))
# #remove it from DF
# N2ad <- N2ad %>% dplyr::select(-shift)
# 
# # get SR HVR calls
# hd_reg <- all_SR_calls %>%
#   dplyr::filter(CHROM==hap_chrom) %>%
#   dplyr::filter((minStart >= hap_start & minStart<= hap_end) | (maxEnd<= hap_end & maxEnd >= hap_start)) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::filter(divSize==max(divSize))
# 
# #get reference-collapsed calls
# hd_collapse <- collapsed_ff %>%
#   dplyr::filter(CHROM==hap_chrom) %>%
#   dplyr::filter((minStart >= hap_start & minStart<= hap_end) | (maxEnd<= hap_end & maxEnd >= hap_start)) %>%
#   dplyr::mutate(minStart=minStart-as.numeric(N2_shift[2]),maxEnd=maxEnd-as.numeric(N2_shift[2])) %>%
#   dplyr::mutate(STRAIN="N2")


# ggplot() +
#   geom_segment(data=N2ad, aes(x=(hdr_start_pos)/1e6,xend=(hdr_end_pos)/1e6, y= sp, yend= sp)) + 
#   geom_rect(data=N2ad,aes(xmin=(newstart + shift)/1e6,xmax=(newend + shift)/1e6,ymin=sp + 0.2, ymax=sp - 0.2,fill=N2))+
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),
#     axis.title = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.text = element_blank()
#   )
# geom_segment(data=N2ad,x=(hdr_start_pos)/1e6,xend=(hdr_end_pos)/1e6, aes(y= sp + 0.1, yend= sp -0.1)) 


########################### PLOT ALL POSSIBLE HAP ##############################

# #gene positions
# plotCoords2 <- rbind(as_tibble(N2ad %>% dplyr::select(-sp)), plotCoords) %>%
#   dplyr::filter(STRAIN=="N2" | STRAIN=="CB4856") %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(minStart=min(newstart,na.rm = T),maxEnd=max(newend,na.rm = T)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(midpoint=(minStart+maxEnd)/2) %>%
#   dplyr::mutate(centeredStart=newstart-midpoint,centeredEnd=newend-midpoint) %>%
#   dplyr::filter(!is.na(seqid)) %>%
#   dplyr::left_join(aliases,by=c("N2"="seqname")) %>%
#   dplyr::mutate(y=ifelse(STRAIN=="N2",2,1))
# 
# #center line positions
# hlines <- plotCoords2 %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(hlineStart=min(centeredStart)-1000,
#                 hlineEnd=max(centeredEnd)+1000) %>%
#   dplyr::select(STRAIN,hlineStart,hlineEnd,y)
# 
# #lab positions
# labs <- plotCoords2 %>%
#   dplyr::mutate(labStart=min(centeredStart)+0.2*(min(centeredStart))) %>%
#   dplyr::select(STRAIN,labStart,seqid,y) %>%
#   dplyr::distinct(seqid,.keep_all = T) 
# 
# regDef_WI <- plotCoords2 %>%
#   #dplyr::filter((!STRAIN=="N2") & (N2=="N2.F19B10.9" | N2=="N2.F40H7.5")) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(regStart=min(centeredStart)+offset,regEnd=max(centeredEnd)-offset) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(STRAIN,.keep_all = T) %>%
#   dplyr::select(STRAIN,regStart,regEnd)
# 
# 
# segments <- plotCoords2 %>%
#   dplyr::filter(STRAIN=="N2") %>%
#   dplyr::mutate(x=((centeredEnd-centeredStart)/2)+centeredStart,y=y-0.2) %>%
#   dplyr::select(N2,x,y,locus_name) %>%
#   dplyr::left_join(plotCoords2 %>%
#                      dplyr::filter(STRAIN=="CB4856") %>%
#                      dplyr::mutate(xend=((centeredEnd-centeredStart)/2)+centeredStart,yend=y+0.2) %>%
#                      dplyr::select(N2,xend,yend,locus_name),by="N2") %>%
#   dplyr::filter(!is.na(xend)) %>%
#   dplyr::distinct() %>%
#   dplyr::group_by(locus_name.x) %>%
#   dplyr::mutate(n1=n()) %>%
#   dplyr::mutate(col1=ifelse(n1>1,"grey","black")) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(locus_name.y) %>%
#   dplyr::mutate(n2=n()) %>%
#   dplyr::mutate(col2=ifelse(n2>1,"grey","black")) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(col=ifelse(col1 == "grey" | col2=="grey","grey","black"))
#   
# 
# cbgenes <- plotCoords2 %>%dplyr::filter(STRAIN=="CB4856")  %>% dplyr::select(locus_name)
# 
# plotCoords3 <- plotCoords2 %>% 
#   dplyr::mutate(N2_color=ifelse(locus_name %in% cbgenes$locus_name,"grey","orange")) %>%
#   dplyr::mutate(N2_color=ifelse(locus_name=="glc-1","green",N2_color))
# #plot
# allhap<- ggplot() +
#   #geom_rect(data=regDef_WI,aes(xmin=regStart,xmax=regEnd,ymin=1-0.5,ymax=1+0.5),fill='lightblue')+
#   geom_segment(data=segments,aes(x=x,y=y,xend=xend,yend=yend,color=col)) +
#   geom_segment(data=hlines,aes(x=hlineStart,xend=hlineEnd,y=y,yend=y),color="black") +
#   geom_rect(data=plotCoords3 %>% dplyr::filter(ortho_status==T),aes(xmin=centeredStart,xmax=centeredEnd,ymin=y-0.2,ymax=y+0.2,fill=N2_color),color="black") +
#   geom_rect(data=plotCoords2 %>% dplyr::filter(ortho_status==F),aes(xmin=centeredStart,xmax=centeredEnd,ymin=y-0.2,ymax=y+0.2,fill="blue"),color='black') +
#   geom_text(data=labs,aes(x=labStart,y=y,label=STRAIN)) +
#  
#   #facet_wrap(~factor(STRAIN, levels=c('N2','NIC2','ECA36','ECA396', 'MY2693', 'EG4725','JU2600','JU310','MY2147','JU1400','NIC526','JU2526','QX1794','CB4856','XZ1516','DL238')),ncol=1) +
#   #facet_wrap(~factor(STRAIN, levels=c('N2','NIC2','ECA36' ,'DL238',"CB4856",'EG4725','MY2147','NIC526','JU2600','JU310','JU1400','XZ1516','JU2526',"QX1794","MY2693","ECA396")),ncol=1) +
#   theme(strip.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         legend.position = 'none',
#         panel.background = element_blank()) +
#   scale_color_manual(values = c('blue'='blue',"black"="grey","grey"="black")) +
#   scale_fill_manual(values=c("orange"="#DB6333","blue"="blue","grey"="grey","lightgrey"="lightgrey","green"="green"))
#   #guides(fill = guide_legend(override.aes = list(size = 0.5))) +
#   #labs(fill='Locus name') 
# allhap


# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/glc1_QTL_CBxN2_CNVPAV.png",allhap, device = 'png',dpi=900,width = 13,height = 1.7,units = 'in')
# save.image(file="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/scripts/Ce_geneAnno-sh/asm_gene_set/glc1_vis_image.Rda")
