#This script processes and standardizes the final set piRNA clusters after liftover
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

# Load libraries
library(dplyr)
library(stringr)

#Load data
#Lifted pics - reference coordinates plus original information

pics_lifted<- read.table("../results/cat_reverted_piCs.txt",header=F)
colnames(pics_lifted)<-c("Chrom_R", "Start_R", "End_R","Strains","piC_number", 
                       "Chrom_RF", "Start_RF", "End_RF", "discard", "Chrom_S", "Start_S", 
                       "End_S","Strain", "infoRPKM","infoReps", "program","discard2" )

#Chrom_R:Start_R-End_R: final reference coordinates of piCs after cross strains intersection 
#Chrom_RF:Start_RF-End_RF: reference coordinates of piCs obtained for each strain
#Chrom_S:Start_S-End_S: strain reference of piCs

pics_lifted<- pics_lifted%>%
  select(piC_number, Chrom_R, Start_R, End_R, Chrom_S, Start_S, End_S,
         Strain, infoRPKM, infoReps, program,
         Chrom_RF, Start_RF, End_RF)

# Removing replicates flagged during quality control
# Replicates removed: SGA12-1, SGA12-2, SGA12-3, SGA34-1, SGA34-2

remove_patterns <- c("SGA12-1", "SGA12-2", "SGA12-3", "SGA34-1", "SGA34-2")
pics_lifted_filtered <- pics_lifted %>%
  rowwise() %>%
  mutate(
    infoRPKM = {
      parts <- str_split(infoRPKM, ";\\s*")[[1]]
      parts <- parts[!str_detect(parts, str_c(remove_patterns, collapse = "|"))]
      if (length(parts) == 0) NA_character_ else str_c(parts, collapse = "; ")
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(infoRPKM))

#cleans and standardizes the program column
pics_lifted_filtered <- pics_lifted_filtered %>%
  rowwise() %>%
  mutate(program = paste(unique(trimws(unlist(strsplit(program, ",")))), collapse = ",")) %>%
  ungroup()
pics_lifted_filtered <- pics_lifted_filtered %>%
  mutate(program = if_else(program == "coverage,protrac", "protrac,coverage", program))

#cleans and standardizes the infoReps column
pics_lifted_filtered <- pics_lifted_filtered %>%
  rowwise() %>%
  mutate(
    infoReps = {
      # Step 1: Split the string by comma, remove whitespace and duplicates
      vals <- unique(trimws(unlist(strsplit(infoReps, ","))))
      
      # Step 2: If shared_C is present, remove exclusive_C
      if ("Replicate-shared_C" %in% vals) {
        vals <- vals[!vals %in% "Replicate-exclusive_C"]
      }
      # Step 3: If shared_P is present, remove exclusive_P
      if ("Replicate-shared_P" %in% vals) {
        vals <- vals[!vals %in% "Replicate-exclusive_P"]
      }
      
      # Step 4: Collapse back into a single comma-separated string
      paste(vals, collapse = ",")
    }
  ) %>%
  ungroup()

#creating Status column
pics_lifted_filtered <- pics_lifted_filtered %>%
  group_by(piC_number) %>%
  mutate(
    unique_sample_count = n_distinct(Strain),
    Status = if_else(unique_sample_count > 1, "Shared", "Exclusive"),
    Strains = paste(sort(unique(Strain)), collapse = ",")
  ) %>%
  ungroup()

# save final result
write.table(pics_filtered, file = "../results/summary_piCs_after_liftover.tsv", row.names=FALSE,sep="\t")
