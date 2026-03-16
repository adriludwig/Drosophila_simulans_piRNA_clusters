########################################################################################
########################### Coverage-based protocol - Step 3 ###########################
########################################################################################

#This script implements Step 3 of the coverage-based piRNA cluster prediction pipeline. 
#It integrates coverage metrics, sequence composition (1U bias), and TE annotations 
#to identify high-confidence piRNA clusters.
# Paths correspond to the internal computing environment and should be adapted if the script is reused.


#The workflow:
# 1 - Merges coverage-based cluster coverage metrics with 1U percentage data.
# 2 - Applies filtering to retain high-confidence candidate piCs (RPKM ≥ 25; 1U bias≥ 60%)
# 3 - Merging overlapping or nearby clusters from different replicates within strains (10kb window)
# 4 - Filters clusters by minimum size (≥2 kb).
# 5- Intersects candidate clusters with TE annotations.
# 6 - Retains only clusters overlapping with TEs and produces the final piC annotation table.

# Load libraries
library(ggplot2)
library(tidyverse)
library(scales)
library(dplyr)
options(scipen = 999)

#Load necessary scripts
source("./scripts/merge_windows.R") # merges clusters detected in different replicates within a 10 kb window.
source("./scripts/merge_with_TE.R") # intersects cluster coordinates with TE annotations.

#Load files:
#potential clusters with RPKM/RPM info
pics_RPKM <- read.table("./tables/all_samples_CPM8_coverage_based.bed", header = TRUE, sep = "\t")
colnames(pics_RPKM) <- c("Chromosome","Start","End","sum_CPM_500bp",
                         "Coverage","Cluster_size","RPKM","RPM","Total_Mapped_Reads","Sample")

#potential clusters with percentage of 1U % info
pics_U1per<-read.table("./tables/all_samples_1U_coverage_based.txt", header = FALSE, sep = "\t") 
colnames(pics_U1per) <- c("U1per", "Sample", "Chromosome", "Start", "End")

#TE annotation
te_bed<- read.table("../TE_Curation_Annotation_Final/TE_annot/all_samples_TE_annot_filtered.bed", header = FALSE, sep = "\t", 
                    stringsAsFactors = FALSE, 
                    quote = "", 
                    fill = TRUE, 
                    comment.char = "")
colnames(te_bed) = c("Chrom", "Start", "End", "TE","NA","Strand","Strain","Size")

#formatting data
pics_RPKM$Sample <- as.character(pics_RPKM$Sample)
pics_RPKM$Chromosome <- as.character(pics_RPKM$Chromosome)
pics_RPKM$Start <- as.numeric(pics_RPKM$Start)
pics_RPKM$End <- as.numeric(pics_RPKM$End)
pics_U1per$Sample <- as.character(pics_U1per$Sample)
pics_U1per$Chromosome <- as.character(pics_U1per$Chromosome)
pics_U1per$Start <- as.numeric(pics_U1per$Start)
pics_U1per$End <- as.numeric(pics_U1per$End)

# Merging dataframes pics_RPKM and pics_U1per
merged_data <- merge(pics_RPKM, pics_U1per, by = c("Sample", "Chromosome", "Start", "End"), all.x = TRUE)
merged_data$Strain <- sapply(strsplit(as.character(merged_data$Sample), "-"), `[`, 1)
merged_data<-merged_data%>%
  mutate(
    Strain = gsub("^MSC\\.01$", "MSC01", Strain),            
    Sample = gsub("\\bMSC\\.01\\b", "MSC01", Sample)         
  )


# Removing replicates flagged during quality control
# Replicates removed: SGA12-1, SGA12-2, SGA12-3, SGA34-1, SGA34-2
merged_data<-merged_data%>%
  filter (!Sample %in% c("SGA12-1","SGA12-2","SGA12-3","SGA34-1","SGA34-2"))

# Filtering by RPKM and U1percentage
#Filter by RPKM 
restrictive_RPKM25<-merged_data%>%
  filter(RPKM >=25)
#Filter by percentage of 1U
restrictive<-restrictive_RPKM25%>%
  filter(U1per >= 60)

# Merging piCs by Strains
# merge piCs from the different replicates within a window of 10kb
merged_df <- merge_windows(restrictive)

# Filtering < 2kb
merged_df<-merged_df%>%
  mutate(Cluster_size=End-Start)

merged_df_size<-merged_df%>%
  filter (Cluster_size>=2000)

#Merging with TE annotation
#intersect potential piCs coordinates and the TE annotation
merged_df_TE <- merge_with_TE(merged_df_size, te_bed)

# Filtering out potential piCs with no TEs
merged_df_TE_filtered <- merged_df_TE %>% filter(!is.na(TE))

# Adding information of the program in the replicates
merged_df_TE_filtered$Replicate_RPKM <- mapply(function(rpkm, strain) {
  # Add "C_" before each replicate entry starting with the strain name
  gsub(paste0("\\b(", strain, "-[^;]+)"), "C_\\1", rpkm)
}, merged_df_TE_filtered$Replicate_RPKM, merged_df_TE_filtered$Strain)

#Editing dataframe and saving file
merged_df_TE_filtered<- merged_df_TE_filtered%>%
  mutate(Status = ifelse(grepl(";", Replicate_RPKM), "Replicate-shared_C", "Replicate-exclusive_C"))%>%
  mutate(program= "coverage" )

data_selected <- merged_df_TE_filtered[, c("Chromosome", "Start", "End", "Strain","Replicate_RPKM", "Status", "program")]

write.table(data_selected, "./tables/piCs_coverage_based.txt", sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
