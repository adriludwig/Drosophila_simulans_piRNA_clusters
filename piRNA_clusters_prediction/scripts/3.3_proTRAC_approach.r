########################################################################################
########################### Coverage-based protocol - Step 3 ###########################
########################################################################################
# This script implements Step 3 of the proTRAC approach to predict piRNA clusters.
# It processes and filters piRNA clusters predicted by proTRAC to generate high-confidence piRNA clusters
# Paths correspond to the internal computing environment and should be adapted if the script is reused.


# The workflow:
# 1 - Imports parsed proTRAC cluster tables generated in Step 2.
# 2 - Merging overlapping or nearby clusters from different replicates within strains (10kb window)
# 3- Intersects candidate clusters with TE annotations.
# 4 - Retains only clusters overlapping with TEs and produces the final piC annotation table.

# 4 - Filters clusters by minimum size ≥ 2 kb).
# 5 - Generates a standardized table of filtered piRNA cluster coordinates and associated metrics.
# 5- Intersects candidate clusters with TE annotations.
# 6 - Retains only clusters overlapping with TEs and produces the final piC annotation table.

# Load libraries
library(ggplot2)
library(tidyverse)
library(scales)
library(dplyr)
options(scipen = 999)

#Load necessary scripts
source("../scripts/merge_windows.R")
source("../scripts/merge_with_TE.R")

#Load files:
protrac <- read.table("../tables/all_samples_protrac_clusters.txt", header = TRUE, sep = "\t")

#TE annotation
te_bed<- read.table("../../TE_Curation_Annotation_Final/TE_annot/all_samples_TE_annot_filtered.bed", header = FALSE, sep = "\t", 
                    stringsAsFactors = FALSE, 
                    quote = "", 
                    fill = TRUE, 
                    comment.char = "")
colnames(te_bed) = c("Chrom", "Start", "End", "TE","NA","Strand","Strain","Size")

# Editing dataframe

protrac<-protrac%>%
  mutate(
    Sample = gsub("\\bMSC\\.01\\b", "MSC01", Sample)) # Partial match in full string

protrac <- select(protrac, -c(Source, Feature, Score, Strand, Frame,piRNA_cluster_no, fraction_piRNA_sized_reads))
protrac<-protrac%>%
  mutate(Strain = sub("-.*", "", Sample))

# Removing replicates flagged during quality control
# Replicates removed: SGA12-1, SGA12-2, SGA12-3, SGA34-1, SGA34-2
protrac<-protrac%>%
  filter (!Sample %in% c("SGA12-1","SGA12-2","SGA12-3","SGA34-1","SGA34-2"))

# Merging piCs by Strains-window of 10 kb
merged_replicates <- merge_windows(protrac)

# Merging with TE annotation
merged_replicates_TEs<- merge_with_TE(merged_replicates, te_bed)

# Filtering out potential piCs without any TEs
merged_replicates_TEs_filtered <- merged_replicates_TEs %>% filter(!is.na(TE))

# Adding information of the program in the replicates
merged_replicates_TEs_filtered$Replicate_RPKM <- mapply(function(rpkm, strain) {
  # Add "C_" before each replicate entry starting with the strain name
  gsub(paste0("\\b(", strain, "-[^;]+)"), "P_\\1", rpkm)
}, merged_replicates_TEs_filtered$Replicate_RPKM, merged_replicates_TEs_filtered$Strain)

# Editing dataframe and saving file
merged_replicates_TEs_filtered<- merged_replicates_TEs_filtered%>%
  mutate(Status = ifelse(grepl(";", Replicate_RPKM), "Replicate-shared_P", "Replicate-exclusive_P"))%>%
  mutate(program= "protrac" )

data_selected <- merged_replicates_TEs_filtered[, c("Chromosome", "Start", "End", "Strain","Replicate_RPKM", "Status", "program")]

write.table(data_selected, "./tables/piCs_protrac.txt", sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
