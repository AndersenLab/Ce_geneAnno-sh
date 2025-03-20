library(plyr) # ALWAYS LOAD BEFORE DPLYR
library(dplyr)
library(ggplot2)
library(readr)
library(IRanges)
library(valr)
library(data.table)



data <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/asssemblies/77HDRs.tsv", col_names = c("CHROM","START","END","strain"))

all_regions <- data %>%
  dplyr::arrange(CHROM,START) %>%
  dplyr::group_split(CHROM)

strain_count <- data %>% dplyr::distinct(strain, .keep_all = T)
print(nrow(strain_count))


getRegFreq <- function(all_regions) {
  all_collapsed <- list()
  for (i in 1:length(all_regions)) {
    temp <- all_regions[[i]]
    k=1
    j=1
    while (k==1) {
      print(paste0("chrom:",i,"/iteration:",j))
      checkIntersect <- temp %>% 
        dplyr::arrange(CHROM,START) %>%
        dplyr::mutate(check=ifelse(lead(START) <= END,T,F)) %>%
        dplyr::mutate(check=ifelse(is.na(check),F,check))
      
      #print(nrow(checkIntersect %>% dplyr::filter(check==T)))
      
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
          dplyr::mutate(newStart=min(START)) %>%
          dplyr::mutate(newEnd=max(END)) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(gid,.keep_all = T)  %>%
          dplyr::mutate(START=newStart,END=newEnd) %>%
          dplyr::select(-newEnd,-newStart)
        
        retain <- temp %>%
          dplyr::filter(check==F & lag(check)==F)
        
        temp <- rbind(collapse,retain) %>%
          dplyr::select(-gid,-check)
        
        j=j+1
      }
    }
    print(head(temp))
    all_collapsed[[i]] <- temp
  }
  return(all_collapsed)
}

HDR_collapse_master <- getRegFreq(all_regions)

all_collapsed <- ldply(HDR_collapse_master, data.frame) %>%
  dplyr::select(-strain) 

colnames(all_collapsed) <- c("chrom","start","end")


vcf <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/vcf/updated_genomes/elegans.merged.1kbCOV.5kbALIGN.annotated.final.vcf") #JU2617 is not present (correctly, it is not an isotype ref so we don't have HDRs for it) 
# and I changed PB306 to it's correct name of ECA259


vfilt <- vcf %>%
  dplyr::filter(INFO == "DEL" | INFO == "INS") %>%
  dplyr::filter(end - POS > 50) %>%
  dplyr::filter(`#CHROM` != "MtDNA") %>% 
  dplyr::rename(CHROM=`#CHROM`)

# print(colnames(vfilt))

pos_bed <- vfilt %>% 
  dplyr::select(CHROM,POS, INFO) %>%
  dplyr::rename(start=POS) %>%
  dplyr::mutate(end=start+1) %>% 
  dplyr::rename(chrom=CHROM)



ggplot() + 
  geom_rect(data = all_collapsed %>% dplyr::rename(CHROM = chrom), aes(xmin = start/1e6, xmax = end/1e6, ymin = 0, ymax = 1.75, fill = "hdr")) + 
  geom_rect(data= vfilt %>% dplyr::filter(INFO == "DEL"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 1, ymax = 1.5, fill = "DEL")) +
  geom_rect(data = vfilt %>% dplyr::filter(INFO == "INS"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 0.25, ymax = 0.75, fill = "INS")) +
  # geom_rect(data = )
  facet_wrap(~CHROM, nrow=6, scales = "free_x") + 
  scale_fill_manual(
    name = "",
    values = c("DEL" = "red", 'INS' = 'blue', 'hdr' = 'gray70'),
    breaks = c("Deletions", "Insertion")
  ) + 
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text("Genome position (Mb)", size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line.x = element_line(),
    legend.position = "none",
    strip.text = element_text(size = 14, color = 'black'),
    panel.background = element_blank(),
    # plot.margin = margin(0, 0, 0, 0)
  )


isec <- valr::bed_intersect(pos_bed,all_collapsed)
isec_del <- isec %>% dplyr::filter(INFO.x == "DEL")
isec_ins <- isec %>% dplyr::filter(INFO.x == "INS")


print(nrow(isec_del)) # 14747
print(nrow(isec_ins)) # 17933
print(nrow(isec)) # 32680

h <- nrow(isec)
print(nrow(vfilt %>% dplyr::filter(INFO == "DEL"))) # 42331
print(nrow(vfilt %>% dplyr::filter(INFO == "INS"))) # 59374
print(nrow(vfilt)) # 101705

percent = (h/total) * 100
print(percent) # 


vfilt_longer <- vfilt %>%
  dplyr::select(-end,-REF,-ALT) %>%
  pivot_longer(
    cols = -c(CHROM, POS, INFO),  # Keep chrom, pos, and info as identifiers
    names_to = "sample",          # New column for sample names
    values_to = "genotype"           # New column for sample values
    ) %>%
  dplyr::mutate(genotype=ifelse(genotype=="./.",0,1)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n_alt=sum(genotype)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_alt)) %>%
  dplyr::mutate(rid=rleid(n_alt)) %>%
  dplyr::filter(genotype > 0)

INDEL_plt <- ggplot(vfilt_longer) + geom_rect(aes(xmin=(POS-500)/1e6,xmax=(POS+500)/1e6,ymin=rid+0.7,ymax=rid-0.7,fill=INFO)) + 
  facet_wrap(~CHROM,nrow=1,scales = 'free_x') +
  scale_fill_manual(values = c("DEL" = "red", "INS" = "blue")) +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  xlab("Genome position (Mb)") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  ylab("77 WI genomes")
INDEL_plt

ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/INDELs_77strains_lowerDPI.jpg", INDEL_plt, dpi = 300, width = 14, height = 10)

