# merge_with_TE.R
# Load required libraries
library(dplyr)
library(GenomicRanges)

# Function to merge dataframe with TE annotation
merge_with_TE <- function(df, te_bed) {
  
  merged_list <- list()  # list to collect merged results
  
  for (strain in unique(df$Strain)) {  # Process each strain separately
    strain_df <- df %>% filter(Strain == strain)
    strain_te <- te_bed %>% filter(Strain == strain)
    
    for (chrom in unique(strain_df$Chromosome)) {  # Process each chromosome separately
      subset_df <- strain_df %>% filter(Chromosome == chrom)
      subset_te <- strain_te %>% filter(Chrom == chrom)
      
      # Only proceed if there are TE annotations for this strain and chromosome
      if (nrow(subset_te) > 0) {
        
        # Convert dataframes to GRanges objects
        df_gr <- GRanges(seqnames = subset_df$Chromosome, 
                         ranges = IRanges(start = as.numeric(subset_df$Start), end = as.numeric(subset_df$End)))
        
        te_gr <- GRanges(seqnames = subset_te$Chrom, 
                         ranges = IRanges(start = as.numeric(subset_te$Start), end = as.numeric(subset_te$End)),
                         TE = subset_te$TE)
        
        # Find overlaps
        hits <- findOverlaps(df_gr, te_gr)
        
        # Initialize TE column
        subset_df$TE <- NA
        
        if (length(hits) > 0) {
          matches <- data.frame(
            df_index = queryHits(hits),
            te_index = subjectHits(hits)
          )
          
          # Fill TE information
          for (i in unique(matches$df_index)) {
            matched_tes <- subset_te$TE[matches$te_index[matches$df_index == i]]
            subset_df$TE[i] <- paste(unique(matched_tes), collapse = "; ")
          }
        }
        
      } else {
        # If no TEs for this chromosome and strain, just add empty TE column
        subset_df$TE <- NA
      }
      
      merged_list[[length(merged_list) + 1]] <- subset_df
    }
  }
  
  final_merged_df <- bind_rows(merged_list)
  
  return(final_merged_df)
}
