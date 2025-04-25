library(plyr) # ALWAYS LOAD BEFORE DPLYR
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ape)
library(tidyr)
library(purrr)
# library(valr)
library(data.table)
# install.packages("nls2")
# library(nls2)


gff <- ape::read.gff("/vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3")
all_gff <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/raw_data/assemblies/elegans/gff/20250423_116_N2_allWS_geneAndmRNA_gff.tsv", col_names = c("type", "start", "end", "strand", "attributes"))

all_gff_formatted <- all_gff %>%
  dplyr::mutate(attributes = gsub("ID=","", attributes)) %>%
  dplyr::mutate(attributes_new = str_extract(attributes, "transcript:[^;]+"), parent = str_extract(attributes, "Parent=gene:[^;]+")
  ) %>%
  dplyr::mutate(attributes_new = ifelse(is.na(attributes_new),attributes,attributes_new)) %>%
  dplyr::mutate(parent = ifelse(is.na(parent),attributes,parent)) %>%
  dplyr::select(-attributes) %>%
  dplyr::mutate(parent = gsub("Parent=gene:","",parent)) %>%
  dplyr::mutate(attributes_new = gsub("transcript:","",attributes_new)) %>%
  dplyr::mutate(attributes_new = sub(";P.*", "", attributes_new), parent = sub(".*Parent=([^;]+);?.*", "\\1", parent)) %>%
  dplyr::mutate(parent = sub(".*gene:(WBGene[0-9]+).*", "\\1", parent), attributes_new = sub(".*gene:(WBGene[0-9]+).*", "\\1", attributes_new)) %>%
  dplyr::mutate(parent = gsub(";","",parent), attributes_new = gsub(";","",attributes_new))

transcript_key <- all_gff_formatted %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::select(attributes_new,parent) %>%
  dplyr::rename(transcript = attributes_new, gene = parent)

# write.table(transcript_key,"/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/transcripts_gene_115WI_N2.tsv", quote = F, row.names = F, col.names = F, sep = '\t')

orthogroups_tran <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/prot_115/OrthoFinder/Results_Apr22/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
# orthogroups <- readr::read_tsv("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/orthology/elegans/prot_78/OrthoFinder/Results_Mar26/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
strCol <- colnames(orthogroups_tran)
strCol_c1 <- gsub(".braker.protein","",strCol)
strCol_c2 <- gsub("_WS283.protein","",strCol_c1)
colnames(orthogroups_tran) <- strCol_c2

# print(nrow(orthogroups_tran)) # 64377

#### Counting number of transcripts per orthogroup for each strain #### 
ortho_count_tran <- orthogroups_tran

strCol_c2_u <- strCol_c2[!strCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]

for (i in 1:length(strCol_c2_u)) {
  print(paste0(i,"out of", length(strCol_c2_u)))
  temp_colname = paste0(strCol_c2_u[i], "_count")

  ortho_count_tran <- ortho_count_tran %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strCol_c2_u[i]),", ") + 1)
}


all_relations_tran <- ortho_count_tran %>%
  dplyr::select(HOG, dplyr::contains("_count"))


#### Transcripts converted to genes #### see script: /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/tran_gene.sh
ortho_genes <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0422_genes.tsv")

stCol <- colnames(ortho_genes)
stCol_c1 <- gsub(".braker.protein","",stCol)
stCol_c2 <- gsub("_WS283.protein","",stCol_c1)
colnames(ortho_genes) <- stCol_c2

ortho_count <- ortho_genes

stCol_c2_u <- stCol_c2[!stCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]

for (i in 1:length(stCol_c2_u)) {
  print(paste0(i,"out of", length(stCol_c2_u)))
  temp_colname = paste0(stCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(stCol_c2_u[i]),", ") + 1)
}

all_relations_duplicated_genes <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))

## all_relations and all_relations_g should be identical
identical(all_relations_tran, all_relations_duplicated_genes) # TRUE

# Remove all duplicate genes in a dataframe cell - see script: /vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/tran_gene.sh
ortho_genes_dd <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/orthofinder/N0_0422_genes_dedup.tsv")

strainCol <- colnames(ortho_genes)
strainCol_c1 <- gsub(".braker.protein","",strainCol)
strainCol_c2 <- gsub("_WS283.protein","",strainCol_c1)
colnames(ortho_genes) <- strainCol_c2

ortho_count <- ortho_genes

strainCol_c2_u <- strainCol_c2[!strainCol_c2 %in% c("OG", "HOG", "Gene Tree Parent Clade")]

for (i in 1:length(strainCol_c2_u)) {
  print(paste0(i,"out of", length(strainCol_c2_u)))
  temp_colname = paste0(strainCol_c2_u[i], "_count")
  
  ortho_count <- ortho_count %>%
    dplyr::mutate(!!sym(temp_colname) := stringr::str_count(!!sym(strainCol_c2_u[i]),", ") + 1)
}

all_relations_duplicated_genes <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count"))





##### Fix for plotting genes and not transcripts
count <- ortho_count %>%
  dplyr::select(HOG, dplyr::contains("_count")) %>%
  dplyr::filter(if_all(-1, ~ . %in% "1" | is.na(.))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u)))


n2_gene <- gff %>%
  dplyr::filter(type == "gene") %>%
  tidyr::separate(attributes, into = c("pre","post"), sep = ";sequence_name=") %>%
  tidyr::separate(pre, into = c("ID","name","L1", "L2"), sep = ";") %>%
  tidyr::separate(post, into = c("seqname","etc"), sep = ";biotype=") %>%
  dplyr::mutate(L2 = ifelse(is.na(L2) & grepl("locus=",L1),L1, NA)) %>%
  dplyr::select(-L1, -name, -etc) %>%
  dplyr::mutate(ID = gsub("ID=Gene:","",ID), L2 = gsub("locus=","",L2)) %>%
  dplyr::mutate(L2 = ifelse(is.na(L2),seqname,L2))

n2_tran <- gff %>%
  dplyr::filter(type == "mRNA") %>%
  tidyr::separate(attributes, into = c("pre","post"), sep = ";Name=") %>%
  tidyr::separate(pre, into = c("ID", "parent"), sep = ";Parent=Gene:") %>%
  dplyr::select(-post) %>%
  dplyr::mutate(ID = gsub("ID=Transcript:","",ID)) %>% 
  dplyr::select(ID,parent) %>%
  dplyr::rename(tran_name = ID) %>%
  dplyr::left_join(n2_gene, by = c("parent" = "ID")) %>%
  dplyr::select(-source, -type, -score, -phase) %>%
  dplyr::rename(locus = L2) %>%
  dplyr::select(seqname,everything())

n2_table <- ortho_genes %>% 
  dplyr::select(HOG,N2) %>%
  dplyr::mutate(N2 = gsub("Transcript_","",N2))

ortho_count_wCoord <- count %>%
  dplyr::left_join(n2_table, by = "HOG") %>%
  dplyr::left_join(n2_tran, by = c("N2" = "tran_name")) %>%
  dplyr::filter(!is.na(seqname)) 


#### Plotting classification based on all orthogroups ####
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
  ylab("Percent of orthogroups") +
  xlab("Frequency") +
  # ggtitle("All orthogroups") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + 
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    legend.position = c(0.85, 0.8),
    plot.title = element_text(size=18, face = 'bold', hjust=0.5),
    legend.text = element_text(size=13, color = 'black'),
    axis.text = element_text(size=12, color = 'black')
  )
gs_allOrtho

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/gene_set_allOrtho_115.png", gs_allOrtho, height = 6, width = 12, dpi = 600)

#### Plotting classification based on 1to1 orthogroups ####
# classification_onetoone <- count %>%
#   dplyr::mutate(
#   class = case_when(
#     freq == 1 ~ "core",
#     freq >= 0.95 & freq < 1 ~ "soft-core",
#     freq > private_freq & freq < 0.95 ~ "accessory",
#     freq == private_freq ~ "private",
#     TRUE ~ "undefined"
#   )
# ) %>%
#   dplyr::count(freq, class) %>%
#   dplyr::mutate(percent = (n / sum(n)) * 100)  
#   
# ggplot(data = classification_onetoone,aes(x = freq * 100, y = percent, fill = class)) + 
#   geom_bar(stat = "identity", color = "black", alpha = 0.5) + 
#   scale_fill_manual(values = c(
#     "core" = "steelblue4",
#     "soft-core" = "steelblue1",
#     "accessory" = "seagreen2",
#     "private" = "seagreen4"
#   ), 
#   limits = c("core", "soft-core", "accessory", "private"),  # Manually ordering legend items
#   guide = guide_legend(title = NULL) 
#   ) +
#   ylab("Percent of orthogroups") +
#   xlab("Frequency (%)") +
#   ggtitle("1-to-1 orthogroups") +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Ensure y-axis is in percentage format
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 16),
#     legend.position = c(0.85, 0.8),
#     plot.title = element_text(size=18, face = 'bold', hjust=0.5),
#     legend.text = element_text(size=13, color = 'black'),
#     axis.text = element_text(size=12, color = 'black')
#   )




#### Plotting 1to1 found in N2 reference ####
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


# Adding resolution if genes are found in a HDR or not
ortho_count_wCoord_HDR <- ortho_count_wCoord %>%
  dplyr::rowwise() %>%
  dplyr::mutate(in_HDR = any(
    seqid == all_collapsed$CHROM &
      start >= all_collapsed$START &
      end <= all_collapsed$END
  ))


ggplot(ortho_count_wCoord_HDR %>% dplyr::filter(!seqid=="MtDNA") %>% dplyr::mutate(in_HDR = factor(in_HDR, levels = c(TRUE, FALSE)))) +
  geom_point(aes(x=start/1e6, y = freq * 100, color=in_HDR)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +  
  geom_smooth(aes(x = start / 1e6, y = freq * 100), method = "loess", se = TRUE, color = "blue") +  
  facet_wrap(~seqid, ncol = 6, nrow = 1, scales = "free_x") + 
  xlab("Genome position (Mb)") + 
  ylab("Gene population frequency (%)") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme(
    legend.position = 'none',
    # legend.text = element_text(size = 13, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    axis.title = element_text(size = 16),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    strip.text = element_text(size = 16, face = "bold", color = "black")  # Larger facet header text
  )

#### Plotting all orthogropus ####

 SOME TYPE OF VISUALIZATION OF GENES FOUND IN N2 VERSUS NOT IN N2??






#### Plotting rarefaction ####
private <- all_relations %>%
  dplyr::mutate(across(2:(ncol(.)), ~ ifelse(. >= 1, 1, .))) %>%
  dplyr::mutate(sum = rowSums(across(-1, ~ ., .names = NULL), na.rm = TRUE)) %>%
  dplyr::mutate(freq = (sum / length(strainCol_c2_u))) %>%
  dplyr::mutate(
    class = case_when(
      freq == 1 ~ "core",
      freq >= 0.95 & freq < 1 ~ "soft-core",
      freq > private_freq & freq < 0.95 ~ "accessory",
      freq == private_freq ~ "private",
      TRUE ~ "undefined"
    )
  ) %>%
  dplyr::filter(class == "private")

strains <- private %>%
  dplyr::select(-HOG,-sum,-freq,-class) %>%
  colnames()

private_ordered <- private %>%
  dplyr::select(all_of(strains)) %>%  # Keep only strain columns
  dplyr::summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "strain", values_to = "count") %>%  # Use a string for `names_to`
  dplyr::arrange(desc(count)) %>%
  dplyr::pull(strain)

private_final <- private %>%
  dplyr::select(HOG, all_of(private_ordered), sum, freq, class)

rarefaction <- data.frame(
  num_strains = integer(),
  num_private_orthogroups = integer()
)

for (i in 2:length(strainCol_c2_u)) {
   temp <- private_final %>% dplyr::select(2:i)
   # print(temp)

   private_count <- sum(apply(temp, 1, function(row) sum(!is.na(row)) == 1))
   rarefaction <- rbind(rarefaction, data.frame(num_strains = i-1, num_private_orthogroups = private_count))
}

rfc <- ggplot(data = rarefaction) +
  geom_point(aes(x=num_strains, y = num_private_orthogroups), size=2.5, color = 'magenta3', alpha=0.8) +
  xlab("Genomes") +
  ylab("Private orthogroups") +
  # coord_cartesian(ylim = c(200,1500)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size=14, color = 'black'),
    axis.text = element_text(size=12, color = 'black'),
    panel.border = element_rect(fill = NA),
  )

rfc

# ggsave("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/plots/rarefaction_115.png", height = 5, width = 9, rfc, dpi = 600)


# Fitting a michaelis-menten model to the data ###
m_model <- nls(num_private_orthogroups ~ Smax * (1 - exp(-k * num_strains)), 
               data = rarefaction,
               start = list(Smax = max(rarefaction$num_private_orthogroups), k = 0.01))

params <- coef(m_model)
Smax <- params["Smax"]
k <- params["k"]

# Compute fitted values & slope (derivative)
rarefaction <- rarefaction %>%
  mutate(predicted = predict(m_model),  
         slope = Smax * k * exp(-k * num_strains))

flat_point <- rarefaction %>%
  filter(slope < 0.5) %>%
  slice(1)  

ggplot(data = rarefaction) +
  geom_line(aes(x = num_strains, y = predicted), color = "blue", linewidth = 1, alpha = 0.5) + 
  geom_point(aes(x = num_strains, y = num_private_orthogroups), size = 3, color = 'seagreen4', alpha = 0.9) +
  geom_vline(xintercept = flat_point$num_strains, linetype = "dashed", color = "red") +  
  annotate("text", x = flat_point$num_strains, y = min(rarefaction$num_private_orthogroups), 
           label = paste("Flattening (slope < 0.5) at", flat_point$num_strains, "genomes"), 
           color = "red", hjust = 1.2) +
  xlab("Number of genomes") +
  ylab("Number of private orthogroups") +
  coord_cartesian(ylim = c(200,1500)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 14, color = 'black'),
    axis.text = element_text(size = 12, color = 'black'),
    panel.border = element_rect(fill = NA)
  )