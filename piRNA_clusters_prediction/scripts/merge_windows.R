# merge_windows_2.R
# Load required libraries
library(dplyr)
library(IRanges)
library(tidyr)

# Function to merge entries within a 10kb window for all strains and replicates
merge_windows <- function(df) {
  df <- df %>% arrange(Strain, Chromosome, Start)  # Sort input
  
  merged_list <- list()
  
  for (strain in unique(df$Strain)) {
    strain_df <- df %>% filter(Strain == strain)
    
    for (chrom in unique(strain_df$Chromosome)) {
      subset_df <- strain_df %>% filter(Chromosome == chrom)
      
      if (nrow(subset_df) > 0) {
        
        # Convert to IRanges
        ranges <- IRanges(start = subset_df$Start, end = subset_df$End)
        
        # Merge ranges within 10kb
        groups <- reduce(ranges, min.gapwidth = 10000)
        overlaps <- findOverlaps(ranges, groups)
        
        # Assign cluster IDs
        subset_df$Cluster <- subjectHits(overlaps)
        
        # Summarize each merged window
        merged_entry <- subset_df %>%
          group_by(Strain, Chromosome, Cluster) %>%
          reframe(
            Start = min(Start),
            End = max(End),
            Replicate_RPKM = paste(Sample, RPKM, sep = ":", collapse = "; ")
          ) %>%
          select(-Cluster)
        
        merged_list[[length(merged_list) + 1]] <- merged_entry
      }
    }
  }
  
  final_merged_df <- bind_rows(merged_list)
  
  return(final_merged_df)
}
