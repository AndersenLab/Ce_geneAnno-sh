library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(data.table)
library(stringr)


# want <- c("N2","JU346","ECA2581","ECA36","ECA1825","ECA1409","ECA1761","ECA1769")
want <- c("N2","ECA1409")

#read all pairwise genome coordinate comparisons
transformed_coords <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/nucmer_runs/115_WI_transformed_coords_FIXED.tsv",col_names = F) 
colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN") 
transformed_coords <- transformed_coords %>%
  dplyr::filter(STRAIN %in% want)

#read concatentated gene models of every genome
gffCat1 <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/braker_runs/merged_gff/all_WI_braker.clean.gff", col_names = F)
colnames(gffCat1) <- c("seqid","source","type","start","end","score","strand","phase","attributes","STRAIN")
gffCat2 <- ape::read.gff("/vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3") %>% dplyr::mutate(STRAIN="N2")
gffCat <- rbind(gffCat1,gffCat2) %>% 
  dplyr::filter(STRAIN %in% want)

#read ortholog relationships among gene models
orthos <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/prot_78/OrthoFinder/Results_Mar20/Orthogroups/Orthogroups.tsv")
strainCol <- colnames(orthos)
strainCol_c1 <- gsub(".braker.protein","",strainCol)
strainCol_c2 <- gsub("_WS283.protein","",strainCol_c1)
colnames(orthos) <- strainCol_c2

# contains glc-1
hdr_chrom = "V"
hdr_start_pos = 16067500 
hdr_end_pos = 16073980

#offset lets you explore adjacent regions
offset = 0
hap_chrom = hdr_chrom
hap_start = hdr_start_pos - offset
hap_end = hdr_end_pos + offset 


#use reference coordinates from g2g alginments to pull the contigs that contain the alt haplotypes for the HDR
hap_coords <- transformed_coords %>%
  dplyr::filter((REF == hap_chrom & hap_start >= S1 & hap_start <= E1 ) | 
                  (REF == hap_chrom & hap_end >= S1 & hap_end <= E1) | 
                  (REF == hap_chrom & S1 >= hap_start & E1 <= hap_end)) %>%
  dplyr::mutate(inv=ifelse(S2>E2,T,F))  %>%
  dplyr::mutate(St2=ifelse(inv==T,E2,S2),Et2=ifelse(inv==T,S2,E2))

# naive visualization of g2g alignments for the target region
# multiple contigs may map to the REF region, we need to filter those!
ggplot(hap_coords) + 
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("REF genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

#keep only the contig with the largest extent of alignment with the REF HDR
tigFilt <- hap_coords %>% 
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(nalign = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,HIFI) %>%
  dplyr::mutate(ntig= n()) %>%
  dplyr::mutate(tigsize=sum(L1)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(tigsize == max(tigsize)) %>%
  dplyr::mutate(rangeDiff=(max(S1,E1)-min(S1,E1))-(hap_end-hap_start)) %>%
  dplyr::ungroup()

#important diagnostic plot
#after the optimal contig is selected, we can see if the HDR boundaries are within the alignment
#for proper visualization, the HDR (grey box in plot) needs to be encompassed by the selected contig
#otherwise some haplotypes will be truncated 
#a step to drop genomes with incomplete coverage of the HDR could be added
#have in mind that the alignments look fragmented because of sequence divergence, but the contig is linear in the genome file
ggplot(tigFilt) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') #+

#keep the set of alignments with the largest span (i.e. removes small distant alignments ("jumps"))
tigFilt2 <- tigFilt %>%
  dplyr::arrange(St2) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDiff=lead(St2)-Et2) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 1.5E5,1,0)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(Et2-St2)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup()

ggplot(tigFilt2) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

trim_spacer = 2e4
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
tigTrim <- tigFilt2 %>%
  dplyr::arrange(STRAIN,S1) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rboundDist=max(S1,E1)-hap_end) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & inv==F,(E2-(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & inv==T,(E2+(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer,(E1-(rboundDist-trim_spacer)),E1)) %>%
  dplyr::mutate(lboundDist=hap_start-min(S1,E1)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & inv==F,(S2+(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & inv==T,(S2-(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer,(S1+(lboundDist-trim_spacer)),S1))

ggplot(tigTrim) +
  geom_rect(xmin=hap_start/1e6,xmax=hap_end/1e6,ymin=-Inf,ymax=Inf,fill="lightgrey")+
  #geom_rect(xmin=16219634/1e6,xmax=16221917/1e6,ymin=-Inf,ymax=Inf,aes(fill="glc-1")) +
  geom_segment(aes(x=S1/1e6,xend=E1/1e6,y=S2/1e6,yend=E2/1e6,color=HIFI)) +
  facet_wrap(~STRAIN,scales = 'free') +
  xlab("N2 genome position (Mb)") +
  ylab("WILD contig position (Mb)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.position = 'none') 

#get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
HV_boundary <- tigTrim %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::select(-leadDiff,-inv) %>%
  dplyr::mutate(inv=ifelse(sum(E2)-sum(S2) >0,F,T)) %>%
  dplyr::mutate(refStart=min(S1,E1),refEnd=max(S1,E1)) %>%
  dplyr::mutate(boundStart=min(S2,E2), boundEnd=max(S2,E2)) %>%
  dplyr::distinct(STRAIN, .keep_all = T) %>%
  dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd,inv) %>%
  dplyr::ungroup() %>%
  dplyr::rename(boundChrom=HIFI)


#filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
wild_genes <- gffCat %>%
  dplyr::filter(type=="gene" & !(STRAIN=="N2")) %>%
  dplyr::mutate(attributes=gsub(";","",attributes)) %>% 
  dplyr::mutate(attributes=gsub("ID=","",attributes)) %>%
  dplyr::select(attributes,seqid,start,end,strand,STRAIN) %>%
  dplyr::rename(Name=attributes)

wild_tran <-  gffCat  %>%
  dplyr::filter(type=="mRNA" & !(STRAIN=="N2")) %>%
  tidyr::separate(attributes,into=c("tranname","Parent"),sep=";Parent=") %>%
  dplyr::mutate(Parent=gsub(";","",Parent)) %>%
  dplyr::mutate(tranname=gsub("ID=","",tranname)) %>%
  dplyr::select(tranname,Parent,STRAIN) %>%
  dplyr::left_join(wild_genes,by=c("Parent"="Name","STRAIN")) %>%
  dplyr::left_join(HV_boundary,by="STRAIN")

#extract the REF genes 
N2Start = min(HV_boundary$refStart)
N2End = max(HV_boundary$refEnd)
N2_genes <- gffCat %>%
  dplyr::filter(type=="gene" & STRAIN=="N2") %>%
  dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter(grepl("biotype=protein_coding",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  #dplyr::mutate(seqname=paste0(seqname)) %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name)) %>% 
  dplyr::select(seqid,start,end,strand,Name,rest2,seqname,STRAIN,refStart,refEnd,seqname) %>%
  dplyr::mutate(rest2=ifelse(grepl("locus",rest2),gsub("locus=","",rest2),seqname)) %>%
  dplyr::rename(alias=rest2)

#extract the REF protein-coding transcripts
N2_tran <- gffCat %>%
  dplyr::filter(type=="mRNA" & STRAIN=="N2") %>%
  tidyr::separate(attributes, into=c("ID","Parent","Name","wormpep","locus"),sep=';') %>%
  dplyr::mutate(ID=gsub("ID=Transcript:","",ID)) %>%
  tidyr::separate(ID,into = c("fosmid","tseqID",'tranum'),sep="\\.",remove = F) %>%
  dplyr::mutate(tseqname=paste0(fosmid,".",tseqID,".",tranum)) %>%
  dplyr::mutate(Parent=gsub("Parent=Gene:","",Parent)) %>%
  dplyr::filter(Parent %in% N2_genes$Name) %>%
  dplyr::select(tseqname,Parent) %>%
  dplyr::rename(tranname=tseqname) %>%
  dplyr::mutate(tranname=paste0("Transcript_",tranname)) %>%
  dplyr::left_join(N2_genes,by=c('Parent'='Name'))


N2_tran_reg <- N2_tran %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | (end >= hap_start & end <= hap_end))  %>%
  dplyr::filter(seqid==hap_chrom)


#get gene list
HV_genelist <- N2_tran_reg$tranname 
#get alt gene names/aliases
aliases <- N2_tran %>% dplyr::select(seqname,tranname,alias)

#minor diagnostic plot to visualize the REF loci captured by the HDR
#this is your REF haplotype
ggplot(N2_tran_reg) + geom_rect(aes(xmin=start,xmax=end,ymin=1,ymax=2))


#filter orthologous groups using REF genes
#this will establish your orthology relationships between REF and WILD haplotypes
all_orthos_unnest <- orthos %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::mutate(na_count = rowSums(is.na(.))) %>%
  dplyr::filter(na_count < length(strainCol_c2) - 2) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end),by=c("N2"="tranname"))  %>%
  dplyr::select(-na_count)

filtOrthos <- orthos %>%
  dplyr::filter(grepl(paste(HV_genelist,collapse="|"),N2)) %>%
  dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
  tidyr::unnest(N2) %>%
  dplyr::mutate(N2=trimws(N2)) %>%
  dplyr::left_join(N2_tran_reg %>% dplyr::select(tranname,seqid,seqname,start,end) %>% dplyr::mutate(og_loc="in_region"),by=c("N2"="tranname")) 

inreg_orthos <- filtOrthos %>% dplyr::filter(!is.na(seqid)) 
outreg_orthos <- filtOrthos %>% dplyr::filter(is.na(seqid)) %>% 
  dplyr::select(-seqid,-seqname,-start,-end,-og_loc) %>%
  dplyr::left_join(N2_tran %>% dplyr::select(tranname,seqid,seqname,start,end,refStart,refEnd) %>% 
                     dplyr::mutate(refChrom=hap_chrom) %>%
                     dplyr::mutate(og_loc="out_region"),by=c("N2"="tranname")) %>%
  dplyr::filter(seqid==refChrom) %>%
  dplyr::mutate(start_dist=abs(refStart-end),end_dist=abs(start-refEnd)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(updown=ifelse(end < refStart,"upstream","downstream")) %>%
  dplyr::mutate(status=ifelse(min_dist_bases <2e4,"out_expand","outside"))

if (nrow(outreg_orthos  %>% dplyr::filter(min_dist_bases < 2e4)) > 0) {
  print("WARNING: There is at least one paralog that is within 10 kb of a gene within your defined boundary in N2. Your boundary will be automatically expanded to include:")
  print(outreg_orthos %>% dplyr::select(seqid,seqname,start,end,N2,min_dist_bases) %>% dplyr::filter(min_dist_bases < 2e4))
}

#generate a lookup table (all_ortho_pairs) which contains all pairwise gene orthologs between REF and WILD
orthoList <- list()
orthoList_bound <- list()
orthoList_raw <- list()
strainCol_iter <- strainCol_c2[!strainCol_c2 %in% c("Orthogroup","N2")]

for (i in 1:length(strainCol_iter)) {
  
  id=strainCol_iter[i]
  raw_tmp <- orthos %>%
    dplyr::select(Orthogroup,strainCol_iter[i],N2) %>%
    dplyr::mutate(str=!!sym(strainCol_iter[i])) %>%
    dplyr::mutate(str = strsplit(as.character(str), ",")) %>%
    tidyr::unnest(str) %>%
    dplyr::mutate(str=trimws(str)) %>%
    dplyr::filter(!is.na(N2)) %>%
    dplyr::select(-strainCol_iter[i]) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    dplyr::mutate(has_any_ortho=T) %>%
    dplyr::left_join(wild_tran,by=c("STRAIN","str"="tranname"))
  #dplyr::select(-N2)
  
  orthoList_raw[[i]] <- raw_tmp
  
  print(paste0("Mapped orthologs for ",i,"/",length(strainCol_iter)," strains."))
  tmp <- rbind(inreg_orthos %>% dplyr::mutate(status="within") %>% dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status),outreg_orthos %>% 
                 dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status)) %>%
    dplyr::select(Orthogroup,N2,strainCol_iter[i],og_loc,status) %>%
    dplyr::rename(tmpSel=strainCol_iter[i]) %>%
    dplyr::mutate(newSel = strsplit(as.character(tmpSel), ",")) %>%
    tidyr::unnest(newSel) %>%
    dplyr::mutate(newSel=trimws(newSel)) %>%
    dplyr::select(-tmpSel) %>%
    dplyr::mutate(STRAIN=strainCol_iter[i]) %>%
    tidyr::separate(newSel,into=c("Name","tnum"),sep="\\.",remove = F) %>%
    dplyr::select(Orthogroup,newSel,Name,STRAIN,N2,-tnum,og_loc,status) %>%
    dplyr::rename(tranname=newSel,Parent=Name) %>%
    dplyr::left_join(wild_tran,by=c("tranname","Parent","STRAIN")) 
  
  orthoList[[i]] <- tmp
  boundg <- tmp %>% 
    dplyr::filter(og_loc=="in_region" | status=="out_expand") %>%
    dplyr::select(-og_loc,-status) %>%
    dplyr::filter(seqid==boundChrom) %>%
    dplyr::mutate(og_loc=ifelse(((start >= boundStart & start <= boundEnd) | (end >= boundStart & end <= boundEnd)),"in_region","out_region")) %>%
    dplyr::mutate(start_dist=ifelse(og_loc=="out_region",abs(boundStart-end),NA),end_dist=ifelse(og_loc=="out_region",abs(start-boundEnd),NA)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(min_dist_bases=min(start_dist,end_dist)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(status=ifelse(og_loc=="in_region","within",ifelse(min_dist_bases < 10000 & !is.na(min_dist_bases),"out_expand","outside")))
  check <- boundg %>% dplyr::filter(status=="out_expand")
  
  if (nrow(check) > 0) {
    print(paste0("WARNING: There is at least one paralog that is within 10 kb of a gene within your derived boundary in ",strainCol_iter[[i]],". Your boundary will be automatically expanded to include:"))
    print(check %>% dplyr::select(seqid,tranname,start,end,strand,STRAIN,min_dist_bases) %>% dplyr::distinct(tranname,.keep_all = T))
    
    sorter <- check %>% dplyr::mutate(updown=ifelse(min_dist_bases==start_dist,"upstream","downstream"))
    upstream <- sorter %>% dplyr::filter(updown=="upstream")
    downstream <- sorter %>% dplyr::filter(updown=="downstream")
    
    if (nrow(upstream) > 0) {
      outer_lim <- max(upstream$end)
      inner <- boundg %>% dplyr::filter(status=="within") %>% dplyr::arrange(start) %>% dplyr::filter(start==min(start))
      inner_lim <- min(inner$start)
      seqid_match <- as.character(unique(inner$seqid))
      extension <- raw_tmp %>%
        dplyr::filter(seqid==seqid_match & start > outer_lim & end <inner_lim) %>%
        dplyr::rename(tranname=str) %>%
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,N2,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
        tidyr::unnest(N2) %>%
        dplyr::mutate(N2=trimws(N2))
      
      boundg_inc <- rbind(extension, boundg %>% dplyr::select(-start_dist,-end_dist,-min_dist_bases))
      orthoList_bound[[i]]  <- boundg_inc %>% dplyr::arrange(start)
    } 
    
    if(nrow(downstream) > 0) {
      outer_lim <- min(downstream$end)
      inner <- boundg %>% dplyr::filter(status=="within") %>% dplyr::arrange(start) %>% dplyr::filter(end==max(end))
      inner_lim <- max(inner$end)
      seqid_match <- as.character(unique(inner$seqid))
      extension <- raw_tmp %>%
        dplyr::filter(seqid==seqid_match & start > inner_lim & end < outer_lim) %>%
        dplyr::rename(tranname=str) %>%
        dplyr::select(Orthogroup,tranname,Parent,STRAIN,N2,everything(),-has_any_ortho) %>%
        dplyr::mutate(og_loc="out_region",status="out_extend") %>%
        dplyr::mutate(N2 = strsplit(as.character(N2), ",")) %>%
        tidyr::unnest(N2) %>%
        dplyr::mutate(N2=trimws(N2))
      
      boundg_inc <- rbind(extension, boundg %>% dplyr::select(-start_dist,-end_dist,-min_dist_bases))
      orthoList_bound[[i]]  <- boundg_inc %>% dplyr::arrange(start)
    } 
    
    
  } else {
    orthoList_bound[[i]] <- boundg %>% 
      dplyr::select(-start_dist,-end_dist,-min_dist_bases) %>% dplyr::arrange(start)
  }
}

all_ortho_pairs  <- ldply(orthoList,data.frame) 
all_ortho_pairs_bound_pre <-ldply(orthoList_bound,data.frame) %>% 
  dplyr::filter(!status=="outside") 

corr_jumps <- all_ortho_pairs_bound_pre %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(leadDist=lead(start)-start) %>%
  dplyr::mutate(leadDist=ifelse(is.na(leadDist),0,leadDist)) %>%
  dplyr::mutate(jump=ifelse(lag(leadDist)>5e4 & lag(status)=="within","JUMP","NOJUMP")) %>%
  dplyr::mutate(jump=ifelse(is.na(jump),"NOJUMP",jump)) %>%
  dplyr::mutate(jumpID=rleid(jump)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,jumpID) %>%
  dplyr::mutate(jgroup_size=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(jgroup_size==max(jgroup_size)) %>%
  dplyr::mutate(keep=T)


all_ortho_pairs_bound <- all_ortho_pairs_bound_pre %>%
  dplyr::arrange(STRAIN,start)# %>%
#dplyr::left_join(corr_jumps %>% dplyr::select(STRAIN,Parent,keep),by=c("STRAIN","Parent")) %>%
#dplyr::filter(!is.na(keep))

all_ortho_pairs_raw <- ldply(orthoList_raw,data.frame) %>% dplyr::select(STRAIN,str,has_any_ortho) %>% dplyr::rename(tranname=str) 

new_boundaries_WI <-  all_ortho_pairs_bound %>%
  dplyr::select(seqid,start,end,STRAIN,tranname,Parent) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(tranname,.keep_all = T) %>%
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(Parent,collapse="-")) %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::ungroup() %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

N2_expand <- rbind(inreg_orthos %>% 
                     dplyr::mutate(status="within"),outreg_orthos %>% 
                     dplyr::select(-refStart,-refEnd,-refChrom,-start_dist,-end_dist,-min_dist_bases,-updown)) %>%
  #dplyr::select(Orthogroup,strainCol_iter[i],N2,seqid,seqname,start,end,og_loc,status)) %>%
  dplyr::filter(!status=="outside")

new_boundaries_N2 <- N2_expand %>%
  dplyr::mutate(minStart=min(start), maxEnd=max(end)) %>%
  dplyr::distinct(N2,.keep_all = T) %>% 
  dplyr::filter(start==minStart | end==maxEnd) %>%
  dplyr::mutate(gene2gene=paste(seqname,collapse="-")) %>%
  dplyr::mutate(STRAIN="N2") %>%
  dplyr::distinct(minStart,.keep_all = T) %>%
  dplyr::select(seqid,minStart,maxEnd,STRAIN,gene2gene)

new_boundaries <- rbind(new_boundaries_WI,new_boundaries_N2) %>%
  dplyr::rename(boundStart=minStart,boundEnd=maxEnd)

#find the bound genes for each strain that are not orthologous
boundGenes <- rbind(wild_tran %>% 
                      dplyr::select(-boundChrom,-boundStart,-boundEnd,-REF,-refStart,-refEnd,-inv) %>% 
                      dplyr::mutate(alias=NA),
                    N2_tran %>% dplyr::select(tranname,seqname,STRAIN,seqid,start,end,strand,alias) %>%
                      dplyr::rename(Parent=seqname)) %>%
  dplyr::left_join(new_boundaries,by=c("STRAIN","seqid")) %>%
  dplyr::filter(!is.na(boundStart)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::filter(start >= boundStart & start <= boundEnd) %>%
  dplyr::ungroup()


N2_ad <- boundGenes %>% 
  dplyr::filter(STRAIN=="N2") %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(tranname %in% all_orthos_unnest$N2,T,F)) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(tranname %in% all_ortho_pairs_bound$N2,T,F)) %>%
  dplyr::arrange(start) %>%
  dplyr::group_by(Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup()

reassess_distal <- N2_ad %>% dplyr::filter(has_any_ortho==T & has_bound_ortho==F) 
distal_ortho <- inreg_orthos %>% 
  dplyr::filter(N2 %in% reassess_distal$tranname) %>% 
  dplyr::mutate(comma_count = stringr::str_count(CB4856, ",")+1) %>%
  dplyr::group_by(CB4856) %>%
  dplyr::mutate(comma_count=sum(comma_count)) %>%
  dplyr::filter(comma_count > 1)

N2_ad_corr <- N2_ad %>%
  dplyr::mutate(has_any_ortho=ifelse(tranname %in% distal_ortho$N2,F,has_any_ortho))

g_count <- length(unique(N2_ad_corr$Parent))

# desired_order <- c("N2","ECA1409","ECA1761","JU346","ECA1769","ECA36","ECA1825","ECA2581")
desired_order <- c("N2","ECA1409")

WI_ad <- boundGenes %>% 
  dplyr::filter(!STRAIN=="N2") %>%
  dplyr::left_join(all_ortho_pairs_raw,by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_any_ortho=ifelse(is.na(has_any_ortho),F,has_any_ortho)) %>%
  dplyr::left_join(all_ortho_pairs_bound %>% dplyr::select(tranname,STRAIN,N2,status) %>% dplyr::filter(N2 %in% N2_ad$tranname),by=c("STRAIN","tranname")) %>%
  dplyr::mutate(tr_has_bound_ortho=ifelse(!is.na(status),T,F)) %>%
  dplyr::select(-status,-has_any_ortho) %>% 
  dplyr::rename(N2_name=N2) %>%
  dplyr::left_join(aliases %>% dplyr::select(-seqname),by=c("N2_name"="tranname")) %>%
  dplyr::mutate(alias.x=alias.y) %>%
  dplyr::select(-alias.y,-N2_name) %>%
  dplyr::rename(alias=alias.x) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(has_any_ortho = any(tr_has_any_ortho)) %>% 
  dplyr::mutate(has_bound_ortho = any(tr_has_bound_ortho)) %>%
  dplyr::select(-tr_has_any_ortho,-tr_has_bound_ortho) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(HV_boundary %>% dplyr::select(STRAIN,inv),by="STRAIN") %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(n_gene=n_distinct(Parent)) %>%
  dplyr::mutate(start_sort = if_else(rep(first(inv), dplyr::n()), -start, start)) %>%
  dplyr::arrange(start_sort, .by_group = TRUE) %>%
  dplyr::mutate(first_gene=first(alias)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(STRAIN = factor(STRAIN, levels = desired_order)) %>%
  dplyr::group_by(STRAIN) %>%
  arrange(first_gene, n_gene, .by_group = TRUE) %>% 
  mutate(order_gene = first(first_gene),
         order_num = first(n_gene)) %>% 
  ungroup() %>%
  # arrange(desc(order_gene), desc(order_num)) %>%
  # arrange(desc(order_gene), order_num) %>%
  dplyr::arrange(desc(STRAIN)) %>%
  select(-order_gene, -order_num) %>%
  dplyr::mutate(g_diff = abs(n_gene-g_count)) %>%
  dplyr::mutate(y_pos=rleid(STRAIN)) %>%
  dplyr::select(-g_diff,-start_sort,-inv,-first_gene,-n_gene) 

N2_ad_corr <- N2_ad_corr %>%
  dplyr::mutate(y_pos=max(WI_ad$y_pos)+1)


all_ad <- rbind(N2_ad_corr,WI_ad) %>% 
  dplyr::arrange(STRAIN,start) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(gen_pos=rleid(Parent)) %>%
  dplyr::mutate(shift=min(start)) %>%
  dplyr::mutate(end=end-min(start),start=start-min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(has_any_ortho==T & has_bound_ortho ==T,2,ifelse(has_any_ortho==T,1,0))) %>%
  dplyr::left_join(HV_boundary %>% dplyr::select(STRAIN,inv),by="STRAIN") %>%
  dplyr::mutate(inv=ifelse(is.na(inv),F,inv)) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(bound_corr=max(end)) %>%
  dplyr::mutate(start=ifelse(inv==T,abs(start-bound_corr),start)) %>%
  dplyr::mutate(end=ifelse(inv==T,abs(end-bound_corr),end))


hlines <- new_boundaries %>% 
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,shift) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN")

segments <- all_ortho_pairs_bound %>%
  dplyr::select(STRAIN,Parent,start,end,N2,strand) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::distinct(Parent,N2,.keep_all = T) %>%
  dplyr::left_join(N2_tran %>% 
                     dplyr::rename(start_N2=start,end_N2=end,strand_N2=strand,chrom_N2=seqid,N2id=STRAIN) %>% 
                     dplyr::select(tranname,chrom_N2,start_N2,end_N2,strand_N2,alias,seqname,N2id),
                   by=c("N2"="tranname")) %>% 
  dplyr::filter(N2 %in% N2_ad$tranname) %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by="STRAIN") %>%
  dplyr::rename(WI_y_pos=y_pos) %>%
  dplyr::left_join(all_ad %>% dplyr::select(STRAIN,y_pos) %>% dplyr::distinct(STRAIN,.keep_all = T),by=c("N2id"="STRAIN")) %>%
  dplyr::rename(N2_y_pos=y_pos) %>%
  dplyr::mutate(N2_shift=min(start_N2),WI_shift=min(start)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(WI_x_pos=(start+((end-start)/2))-WI_shift, N2_x_pos=(start_N2+((end_N2-start_N2)/2)-N2_shift)) %>%
  dplyr::mutate(WI_y_pos=WI_y_pos+0.2,N2_y_pos=N2_y_pos-0.2) %>%
  dplyr::distinct(STRAIN,Parent,seqname,.keep_all = T) %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::mutate(n1=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN,seqname) %>%
  dplyr::mutate(n2=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(col=ifelse(n1>1 | n2>1,"multi_copy","single_copy")) 

plot_ad <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) #%>%


all_hap <- ggplot() +
  geom_segment(data=hlines,aes(x=boundStart-shift,xend=boundEnd-shift,y=y_pos,yend=y_pos))+
  geom_rect(data=plot_ad, aes(xmin=start,xmax=end,ymin=y_pos+0.2,ymax=y_pos-0.2,fill=class),color="black") +
  scale_y_continuous(breaks = hlines$y_pos, labels = hlines$STRAIN) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values=c("has_local_ortho"="grey","has_distal_ortho"="black","no_known_ortho"="red","no_known_allelic_CB"="blue")) +
  scale_x_continuous(expand = c(0.01,0)) +
  ylab("77 isotype strains")
all_hap

plot_ad <- all_ad %>%
  dplyr::group_by(STRAIN,Parent) %>%
  dplyr::filter(col==max(col)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(STRAIN,Parent,.keep_all = T) %>%
  dplyr::mutate(class=ifelse(col==0,"no_known_ortho",ifelse(col==1,"has_distal_ortho","has_local_ortho"))) %>%
  dplyr::mutate(alise = as.character(alias))

# Exclude "non-ortho" from the trapezium joining
plot_ad_filtered <- plot_ad %>% 
  dplyr::mutate(alias=ifelse(is.na(alias),"non-ortho",alias)) %>%
  filter(alias != "non-ortho")

# Join filtered data frames for many-to-many connections
trapeziums <- inner_join(
  plot_ad_filtered, plot_ad_filtered,
  by = "alias",
  suffix = c("_upper", "_lower"),
  relationship = "many-to-many"
) %>% 
  filter(y_pos_upper - y_pos_lower == 1)

# Create trapezium polygons using min/max for x-coordinates so that start/end orientation is corrected.
trapezium_polys <- trapeziums %>% 
  rowwise() %>%
  do({
    # Calculate corrected x coordinates for the upper rectangle
    x_left_upper <- min(.$start_upper, .$end_upper)
    x_right_upper <- max(.$start_upper, .$end_upper)
    
    # Calculate corrected x coordinates for the lower rectangle
    x_left_lower <- min(.$start_lower, .$end_lower)
    x_right_lower <- max(.$start_lower, .$end_lower)
    
    data.frame(
      alias = .$alias,
      group = paste(.$alias, .$y_pos_upper, sep = "_"),
      x = c(x_left_upper, x_right_upper, x_right_lower, x_left_lower),
      y = c(.$y_pos_upper - 0.2,  # bottom edge of the upper rectangle
            .$y_pos_upper - 0.2,
            .$y_pos_lower + 0.2,  # top edge of the lower rectangle
            .$y_pos_lower + 0.2)
    )
  }) %>%
  ungroup()

# Extract unique aliases at y_pos 77 in order of increasing start position

ordered_aliases <- plot_ad %>%
  # filter(y_pos == 77) %>%
  filter(y_pos == max(WI_ad$y_pos)) %>%
  arrange(start) %>%
  pull(alias) %>%
  unique()

# Reorder the factor levels so that the legend follows the ordered aliases
plot_ad <- plot_ad %>%
  mutate(alias = factor(alias, levels = ordered_aliases))

# Also update any other data frames with alias info, e.g. trapezium_polys:
trapezium_polys <- trapezium_polys %>%
  mutate(alias = factor(alias, levels = ordered_aliases))

# Shuffle the assignment of colors to the ordered aliases
set.seed(123)  # for reproducibility
shuffled_aliases <- sample(ordered_aliases)

# Generate colors using hcl.colors() for the shuffled aliases
default_colors <- setNames(hcl.colors(length(shuffled_aliases), "Dark 3"), shuffled_aliases)

# But to keep the legend order as ordered_aliases, we re-map these colors back:
final_colors <- default_colors[ordered_aliases]

# Optionally, if you have the "non-ortho" alias (or any other), add it explicitly:
final_colors <- c(final_colors, "non-ortho" = "black")

# Then apply final_colors in your ggplot scale:
all_hap_bg <- ggplot() +
  geom_segment(data = hlines, 
               aes(x = boundStart - shift, xend = boundEnd - shift, y = y_pos, yend = y_pos)) +
  geom_polygon(data = trapezium_polys, 
               aes(x = x, y = y, group = group, fill = alias)) +
  geom_rect(data = plot_ad, 
            aes(xmin = start, xmax = end, ymin = y_pos + 0.2, ymax = y_pos - 0.2, fill = alias),color="black") +
  theme(
    legend.title = element_blank(),
    # legend.position = 'none',
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, color = 'black', face = 'bold'),
    axis.ticks = element_blank()
  ) +
  scale_fill_manual(values = final_colors) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(breaks = hlines$y_pos, labels = hlines$STRAIN, expand = c(0, 0)) +
  ylab("77 isotype strains")

all_hap_bg

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/haplotype_plotter.png", dpi = 600, all_hap_bg, width = 15, height = 8)
ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/haplotype_plotter_threeGenes.png", dpi = 600, all_hap_bg, width = 15, height = 8)

