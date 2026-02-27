library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ape)
library(valr)
library(stringr)
library(data.table)
library(cowplot)


# Testing out Nic's code!
SOI = "ECA3088"

nuc <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/synteny_vis/elegans/nucmer_aln_WSs/142_nucmer_ECA741CGC1.tsv", col_names = c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")) %>%
  dplyr::filter(STRAIN != "ECA396") %>%
  dplyr::filter(L2 > 3e3) %>%
  dplyr::group_by(STRAIN,REF) %>%
  dplyr::arrange(STRAIN,REF,S1) %>%
  dplyr::mutate(leadS1=lead(S1),leadE1=lead(E1),leadS2=lead(S2),leadE2=lead(E2),lagS1=lag(S1),lagE1=lag(E1),lagS2=lag(S2),lagE2=lag(E2)) %>%
  dplyr::ungroup() 

hdrs <- readr::read_tsv("/vast/eande106/data/c_elegans/WI/divergent_regions/20250625/20250625_c_elegans_divergent_regions_strain.bed", col_names = c("chrom", "start", "end", "strain"))

strain_hdrs <- hdrs %>%
  dplyr::arrange(chrom,start) %>%
  dplyr::filter(strain%in%unique(nuc$STRAIN)) %>%
  #dplyr::filter(strain == SOI
  dplyr::rename(STRAIN = strain, CHROM = chrom) %>% 
  dplyr::mutate(REF="N2")


calls <- as.data.table(strain_hdrs)[, rowid := .I]
calls <- calls[, .(rowid, region_source = REF, genome = STRAIN, contig = CHROM,
                   start, end)]

nuc_dt <- as.data.table(nuc)[, idx := .I]
nuc_dt <- nuc_dt[, .(idx, genome = STRAIN, contig = REF,
                     start = S1, end = E1,
                     orig_start=S2,orig_end=E2,
                     L1, L2,
                     leadS1,leadE1,leadS2,leadE2,
                     lagS1,lagE1,lagS2,lagE2,
                     refchrom=HIFI,refstart=pmin(S2, E2),refend=pmax(S2, E2))]

# Remove old values in case you're re-running in the same session
#nuc_dt[, `:=`(region_id = NULL, region_source = NULL)]
calls[, `:=`(start = as.integer(start), end = as.integer(end))]
nuc_dt[, `:=`(start = as.integer(start), end = as.integer(end))]
# Set proper keys
setkey(calls, genome, contig, start, end)
setkey(nuc_dt, genome, contig, start, end)

# Rerun safe and clean overlap
matched <- foverlaps(
  x = calls,
  y = nuc_dt,
  type = "any",
  nomatch = 0
)

matched_df <- as.data.frame(matched) %>%
  dplyr::mutate(hdr_chrom=contig) %>%
  dplyr::mutate(INV=ifelse(orig_start==refstart,F,T)) %>%
  dplyr::select(refstart,refend,orig_start,orig_end,L1,L2,refchrom,contig,genome,INV,start,end,rowid,region_source,hdr_chrom,i.start,i.end,leadS1,leadE1,leadS2,leadE2,lagS1,lagE1,lagS2,lagE2) %>%
  dplyr::rename(S1=start,E1=end,S2=orig_start,E2=orig_end,REF=refchrom,HIFI=contig,HIFI_strain=genome,St2=refstart,Et2=refend,group_id=rowid,hdr_strain=region_source,hdr_start=i.start,hdr_end=i.end) %>% 
  dplyr::arrange(group_id,S1) %>%
  dplyr::mutate(HDRid = paste0(hdr_strain,hdr_chrom,hdr_start,hdr_end))



ggplot(matched_df %>% dplyr::filter(HDRid=="N2I650000660000")) + 
  geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
  geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2)) +
  xlab("REF")+
  ylab("WI")

tigFilt2 <- matched_df %>%
  dplyr::arrange(group_id,S1) %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(leadDiff=lead(S1)-E1) %>%
  dplyr::mutate(jump=ifelse(leadDiff > 5E4,1,0)) %>%
  dplyr::mutate(leadDiff=ifelse(is.na(leadDiff),0,leadDiff)) %>%
  dplyr::mutate(run_id = cumsum(c(1, head(jump, -1)))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group_id,run_id) %>%
  dplyr::mutate(gsize=n()) %>%
  dplyr::mutate(len=abs(E1-S1)) %>%
  dplyr::mutate(sumlen=sum(len)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group_id) %>%
  dplyr::filter(sumlen==max(sumlen)) %>%
  dplyr::select(-gsize) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group_id,S1) 


trim_spacer = 1e3
#trims long alignments to the focal region (i.e. hap_start to hap_end, but transformed to the other genome)
tigTrim <- tigFilt2 %>%
  dplyr::group_by(group_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(scale_distortion = ((L2 - L1)/L1)) %>%
  dplyr::mutate(lboundDist=hdr_start-min(S1,E1)) %>%
  dplyr::mutate(rboundDist=max(S1,E1)-hdr_end) %>%
  #dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==F,(E1-(rboundDist-trim_spacer)),E1)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==F ,(E1-(rboundDist-trim_spacer+(rboundDist*scale_distortion))),E1)) %>%
  dplyr::mutate(E1=ifelse(rboundDist>trim_spacer & INV==T ,(E1-(rboundDist-trim_spacer+(rboundDist*scale_distortion))),E1)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & INV==F,(E2-(rboundDist-trim_spacer)),E2)) %>%
  dplyr::mutate(E2=ifelse(rboundDist>trim_spacer & INV==T,(E2+(rboundDist-trim_spacer)),E2)) %>%

  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer & INV==F ,(S1+(lboundDist-trim_spacer+(lboundDist*scale_distortion))),S1)) %>%
  
  dplyr::mutate(S1=ifelse(lboundDist>trim_spacer & INV==T ,(S1+(lboundDist-trim_spacer+(lboundDist*scale_distortion))),S1)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & INV==F,(S2+(lboundDist-trim_spacer)),S2)) %>%
  dplyr::mutate(S2=ifelse(lboundDist>trim_spacer & INV==T,(S2-(lboundDist-trim_spacer)),S2)) %>%
  dplyr::ungroup()

cowplot::plot_grid(
ggplot(tigFilt2%>% dplyr::filter(HDRid=="N2I650000660000")) + 
  geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
  geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2,color=INV)) +
  xlab("REF")+
  ylab("WI"),

ggplot(tigTrim%>% dplyr::filter(HDRid=="N2I650000660000")) + 
  geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
  geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2)) +
  xlab("REF")+
  ylab("WI"))


tigMarkExtend <- tigTrim %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(E2_extend=ifelse(INV==F & E2 == max(E2) & E1 < hdr_end, T, F),
                S2_extend=ifelse(INV==F & S2 == min(S2) & S1 > hdr_start, T,F),
                iE2_extend=ifelse(INV==T & E2 == min(E2) & E1 < hdr_end,T,F),
                iS2_extend=ifelse(INV==T & S2 == max(S2) & S1 > hdr_start, T,F)) %>%
  dplyr::mutate(any_extend=ifelse(E2_extend == T | S2_extend == T | iE2_extend==T | iS2_extend ==T,T,F)) %>%
  dplyr::ungroup()


cowplot::plot_grid(
  ggplot(tigTrim%>% dplyr::filter(HDRid=="N2I12130001218000")) + 
    geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
    geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2,color=INV)) +
    xlab("REF")+
    ylab("WI"),
  
  ggplot(tigMarkExtend%>% dplyr::filter(HDRid=="N2I12130001218000")) + 
    geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
    geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2,color=any_extend)) +
    xlab("REF")+
    ylab("WI"),
  ggplot(tigMarkExtend%>% dplyr::filter(HDRid=="N2I12130001218000")) + 
    geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
    geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2,color=iS2_extend)) +
    xlab("REF")+
    ylab("WI"),
ggplot(tigMarkExtend%>% dplyr::filter(HDRid=="N2I12130001218000")) + 
  geom_rect(aes(xmin=hdr_start,xmax=hdr_end,ymin=Inf,ymax=-Inf)) + 
  geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2,color=iE2_extend)) +
  xlab("REF")+
  ylab("WI"),nrow=1)

tigToExtend <- tigMarkExtend %>% 
  dplyr::filter(any_extend==T) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(extend_length_WI_lead=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end, min(leadS2,leadE2)-E2,
                                             ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end, min(leadS2,leadE2)-S2,NA)),
                extend_length_WI_lag=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start, S2-max(lagS2,lagE2),
                                            ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start, E2-max(lagS2,lagE2),NA)),
                extend_length_REF_lead=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end, 
                                              ifelse(leadS1 >= E1,leadS1-E1,ifelse(leadE1>=E1,leadE1-E1,NA)),
                                              ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end, ifelse(S1>=leadE1,S1-leadE1,ifelse(S1>=leadS1,S1-leadS1,NA)),NA)),
                extend_length_REF_lag=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start, 
                                             ifelse(lagE1<=S1,S1-lagE1,ifelse(lagS1<=S1,S1-lagS1,NA)),
                                             ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start, ifelse(lagS1>=E1,lagS1-E1,ifelse(lagE1>=E1,lagE1-E1,NA)),NA))) %>%
  dplyr::ungroup()


extendDat <- rbind(tigToExtend %>% 
                     dplyr::select(extend_length_REF_lead,extend_length_WI_lead) %>% 
                     dplyr::rename(extend_length_WI=extend_length_WI_lead,extend_length_REF=extend_length_REF_lead),
                   tigToExtend %>% 
                     dplyr::select(extend_length_REF_lag,extend_length_WI_lag) %>% 
                     dplyr::rename(extend_length_WI=extend_length_WI_lag,extend_length_REF=extend_length_REF_lag)) %>%
  dplyr::filter(!is.na(extend_length_WI) & !is.na(extend_length_REF))

sc <- ggplot(data=extendDat) + 
  geom_point(aes(x=extend_length_REF/1e3,y=extend_length_WI/1e3),size=1) + 
  geom_rect(xmin=0,xmax=100,ymin=-Inf,ymax=Inf,fill=NA,color="grey",linetype="dashed")+
  geom_rect(xmin=Inf,xmax=-Inf,ymin=0,ymax=100,fill=NA,color="grey",linetype="dashed")+
  theme_classic() + 
  xlab("REF extension distances (kb)") + 
  ylab("WI extension distances (kb)") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_x_continuous(expand = c(0.01,0))

h1 <- ggplot(data=extendDat) + 
  geom_histogram(aes(x=extend_length_REF/1e3),binwidth = 1) + 
  theme_classic() + 
  xlab("") + 
  ylab("count") + 
  coord_cartesian(xlim=c(0,100)) +
  scale_y_continuous(expand = c(0.01,0),
                     labels = function(y) y / 1000,
                     name = "count (thousand)") +
  scale_x_continuous(expand = c(0.01,0))

h2 <- ggplot(data=extendDat) + 
  geom_histogram(aes(y=extend_length_WI/1e3),binwidth = 1) + 
  theme_classic() + ylab("") + 
  coord_cartesian(ylim=c(0,100)) +
  scale_x_continuous(labels = function(x) x / 1000,
                     expand = c(0.01, 0),
                     name = "count (thousand)") +
  scale_y_continuous(expand = c(0.01,0))


middle_row <- cowplot::plot_grid(
  sc,       
  h2, 
  ncol = 2,
  rel_widths = c(4, 1), 
  align = "hv"
)

empty_plot <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void() 

top_row <- plot_grid(
  empty_plot,
  h1,
  empty_plot,
  ncol = 3,
  rel_widths = c(0.09,4, 1),
  align = "hv"
)

final_plot <- cowplot::plot_grid(
  top_row,      
  middle_row,    
  ncol = 1,
  rel_heights = c(1, 4),  
  align = "v"
)
final_plot

tigExtensions <- rbind(tigToExtend,tigMarkExtend %>% dplyr::filter(any_extend==F) %>% dplyr::mutate(extend_length_WI_lead=NA,extend_length_REF_lead=NA,extend_length_WI_lag=NA,extend_length_REF_lag=NA))




tigExtended_50kb <- tigExtensions %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(E1=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, ifelse(leadS1 >= E1,leadS1,ifelse(leadE1>=E1,leadE1,E1)),E1)) %>%
  dplyr::mutate(E1=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, ifelse(lagS1>=E1,lagS1,ifelse(lagE1>=E1,lagE1,E1)),E1)) %>%
  dplyr::mutate(S1=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, ifelse(lagE1<=S1,lagE1, ifelse(lagS1<=S1,lagS1,S1)),S1)) %>%
  dplyr::mutate(S1=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, ifelse(S1>=leadE1,leadE1,ifelse(S1>=leadS1,leadS1,S1)),S1)) %>%
  dplyr::mutate(E2=ifelse(E2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, min(leadS2,leadE2),E2)) %>%
  dplyr::mutate(E2=ifelse(iE2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, max(lagS2,lagE2),E2)) %>%
  dplyr::mutate(S2=ifelse(S2_extend==T & max(lagS2,lagE2) < hdr_start & !is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag) & extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4, max(lagS2,lagE2),S2)) %>%
  dplyr::mutate(S2=ifelse(iS2_extend==T & min(leadS2,leadE2) > hdr_end & !is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead) & extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4, min(leadS2,leadE2),S2)) %>% 
  dplyr::ungroup()

id="N2I12130001218000"

ggplot(rbind(tigExtensions %>% dplyr::filter(HDRid == id) %>% dplyr::mutate(type="original"),tigExtended_50kb %>% dplyr::filter(HDRid == id) %>% dplyr::mutate(type="extended_max50kb"))) +
  geom_rect(aes(ymin=-Inf,ymax=Inf,xmin=hdr_start,xmax=hdr_end),fill="lightgrey") +
  geom_segment(aes(x=leadS1,xend=leadE1,y=leadS2,yend=leadE2, color="outreg_lead")) +
  geom_segment(aes(x=lagS1,xend=lagE1,y=lagS2,yend=lagE2, color="outreg_lag")) +
  geom_segment(aes(x=S1,xend=E1,y=S2,yend=E2, color="inreg_trimmed")) +
  xlab("REF coords") +
  ylab("WI coords") +
  facet_wrap(~as.factor(type),scales = 'free')

# counts25kb <- tigExtended_25kb %>%
#   dplyr::filter(
#     any_extend == TRUE,
#     (!is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag)) | (!is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead)),
#     (extend_length_WI_lag < 2.5e4 & extend_length_REF_lag < 2.5e4) | (extend_length_WI_lead < 2.5e4 & extend_length_REF_lead < 2.5e4)
#   ) %>%
#   dplyr::group_by(hdr_strain, HDRid) %>%
#   dplyr::summarise(has_extension = any(any_extend), .groups = "drop") %>%  # optional, since filter already ensures this
#   dplyr::group_by(hdr_strain) %>%
#   dplyr::summarise(count_true = n(), .groups = "drop") %>%
#   dplyr::mutate(window_size = "50kb")

counts50kb <- tigExtended_50kb %>%
  dplyr::filter(
    any_extend == TRUE,
    (!is.na(extend_length_WI_lag) & !is.na(extend_length_REF_lag)) | (!is.na(extend_length_WI_lead) & !is.na(extend_length_REF_lead)),
    (extend_length_WI_lag < 5e4 & extend_length_REF_lag < 5e4) | (extend_length_WI_lead < 5e4 & extend_length_REF_lead < 5e4)
  ) %>%
  dplyr::group_by(hdr_strain, HDRid) %>%
  dplyr::summarise(has_extension = any(any_extend), .groups = "drop") %>%  # optional, since filter already ensures this
  dplyr::group_by(hdr_strain) %>%
  dplyr::summarise(count_true = n(), .groups = "drop") %>%
  dplyr::mutate(window_size = "50kb")

hdr_counts <- tigTrim %>%
  dplyr::group_by(hdr_strain) %>%
  dplyr::summarise(num_unique_HDRid = n_distinct(HDRid), .groups = "drop")

hdr_transformed_orig <- tigTrim %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(
    S2 = min(S2, E2),
    E2 = max(S2, E2),
    across(
      .cols = -c(S1, E1, S2, E2, St2, Et2),
      .fns = dplyr::first
    ),
    .groups = "drop"
  )

hdr_transformed_50ext <- tigExtended_50kb %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(
    S2 = min(S2, E2),
    E2 = max(S2, E2),
    across(
      .cols = -c(S1, E1, S2, E2, St2, Et2),
      .fns = dplyr::first
    ),
    .groups = "drop"
  )


test<-hdr_transformed_50ext %>% 
  dplyr::mutate(mode="TRANS") %>%
  dplyr::select(REF,S2,E2,hdr_strain,HIFI_strain,HIFI,hdr_start,hdr_end) %>%
  dplyr::rename(HIFI=REF,minStart=S2,maxEnd=E2,REF=hdr_strain,STRAIN=HIFI_strain,CHROM=HIFI) %>% 
  dplyr::mutate(divSize=abs(maxEnd-minStart)) %>%
  dplyr::mutate(divSizeRef=abs(hdr_start-hdr_end)) %>%
  dplyr::mutate(divDiff=divSize-divSizeRef) %>%
  dplyr::filter(divSize >= 5e3) %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ncalls=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(REF) %>%
  dplyr::arrange(REF,desc(ncalls),STRAIN,CHROM,minStart) %>%
  dplyr::mutate(sorter=paste0(ncalls,STRAIN)) %>%
  dplyr::mutate(rleID=data.table::rleid(sorter)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(STRAIN) %>%
  dplyr::mutate(ystrain=cur_group_id()) %>%
  dplyr::ungroup() 

ggplot(test %>% dplyr::filter(divDiff<5e4 & STRAIN==SOI)) + 
  geom_rect(aes(ymin=minStart/1e6,ymax=maxEnd/1e6,xmin=Inf,xmax=-Inf)) + 
  geom_rect(aes(xmin=hdr_start/1e6,xmax=hdr_end/1e6,ymin=Inf,ymax=-Inf)) + 
  facet_grid(HIFI~CHROM, scales = 'free_x') + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.text.y=element_blank())+
  #strip.text.x = element_blank())  +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("N2 position (Mb)") +
  ylab("WI position (Mb)")



































