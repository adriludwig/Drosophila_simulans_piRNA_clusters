#Step 4 - Merging piCs and 
#This script processes and standardizes the final set piRNA clusters
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

# Load libraries
library(dplyr)
library(stringr)

# Load data
pics <- read.table("./results/pic_outside_gene_2kb_final_all.bed", header = FALSE, sep = "\t") 
colnames(pics) <- c("Chromosome","Start","End","Strain", "Replicate_cov","Replicate_status", "program")

# Removing replicates flagged during quality control
# Replicates removed: SGA12-1, SGA12-2, SGA12-3, SGA34-1, SGA34-2
remove_patterns <- c("SGA12-1", "SGA12-2", "SGA12-3", "SGA34-1", "SGA34-2")
rx_remove <- paste0("^(C_|P_)(", paste(remove_patterns, collapse = "|"), ")(?=:)")

pics_filtered <- pics %>%
  mutate(
    Replicate_cov = str_split(Replicate_cov, "\\s*[;,]\\s*") |>
      lapply(function(x) {
        x <- x[x != ""]
        x <- x[!str_detect(x, rx_remove)]
        if (length(x) > 0) str_c(x, collapse = "; ") else NA_character_
      }) |>
      unlist()
  ) %>%
  filter(!is.na(Replicate_cov))

#Editing dataframe

# Get cluster size
pics_filtered <-pics_filtered%>%
  mutate(Cluster_size=End-Start)

#cleans and standardizes the Replicate_status column
pics_filtered <- pics_filtered %>%
  rowwise() %>%
  mutate(
    Replicate_status = {
      # Step 1: Split the string by comma, remove whitespace and duplicates
      vals <- unique(trimws(unlist(strsplit(Replicate_status, ","))))
      
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

#cleans and standardizes the program column
pics_filtered <- pics_filtered %>%
  rowwise() %>%
  mutate(program = paste(unique(trimws(unlist(strsplit(program, ",")))), collapse = ",")) %>%
  ungroup()
pics_filtered <- pics_filtered %>%
  mutate(program = if_else(program == "coverage,protrac", "protrac,coverage", program))

# save final result
write.table(pics_filtered, file = "./results/summary_piCs_before_liftover.tsv", row.names=FALSE,sep="\t")

