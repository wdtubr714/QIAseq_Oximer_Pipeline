##VALIDATED TARGETS

# Find validated targets and their 3' UTR sequences using a list of differentially expressed miRNAs

# Load the libraries
library(biomaRt)
library(multiMiR)
library(dplyr)
library(readr)
library(tidyr)

# Read file and process
miRNA_file <- "/path/to/file/.csv"
miRNA_list <- read.csv(miRNA_file)

# Retrieve mirnas
mirna_names <- miRNA_list$miRNA

# Remove duplicates
unique_miRNA_names <- unique(mirna_names)

# Create a new data frame with the unique miRNA names
unique_miRNA_df <- data.frame(miRNA = unique_miRNA_names)

# Strictly miRNA column
mirna_names_shortened <- unique_miRNA_df$miRNA

# Plug miRNA's into multiMiR and get validated targets
multimir_results <- get_multimir(org     = 'mmu',
                                 mirna   = mirna_names_shortened,
                                 table   = 'validated',
                                 summary = TRUE)

mirna_targets = as.data.frame(multimir_results@data)

# Function to connect to Ensembl with fallback to different mirrors
connect_to_ensembl <- function() {
  mirrors <- c("www", "useast", "uswest", "asia")
  for (mirror in mirrors) {
    ensembl <- tryCatch({
      useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = mirror)
    }, error = function(e) {
      message("Ensembl site unresponsive, trying ", mirror, " mirror")
      return(NULL)
    })
    if (!is.null(ensembl)) {
      message("Connected to Ensembl using ", mirror, " mirror")
      return(ensembl)
    }
  }
  stop("Unable to connect to any Ensembl mirror")
}

# Run batch
ensembl <- connect_to_ensembl()

# Filter out rows with missing or empty Ensembl IDs
mirna_targets <- mirna_targets[mirna_targets$target_ensembl != "", ]
mirna_targets <- mirna_targets[!is.na(mirna_targets$target_ensembl), ]

# Check the cleaned Ensembl IDs
print("Cleaned Target Ensembl IDs:")
print(unique(mirna_targets$target_ensembl))

# List of cleaned target Ensembl IDs
target_ensembl_ids <- unique(mirna_targets$target_ensembl)

# Define a function to retrieve sequences in batches
get_sequences_in_batches <- function(ids, batch_size = 50, max_retries = 3) {
  total_ids <- length(ids)
  sequence_info <- data.frame()
  
  for (start in seq(1, total_ids, by = batch_size)) {
    end <- min(start + batch_size - 1, total_ids)
    batch_ids <- ids[start:end]
    
    # Retry mechanism
    attempt <- 1
    success <- FALSE
    
    while (attempt <= max_retries && !success) {
      try({
        # Retrieve sequence information for the batch
        batch_sequence_info <- getSequence(id = batch_ids, type = "ensembl_gene_id", seqType = "3utr", mart = ensembl)
        
        # Combine with the main data frame
        sequence_info <- rbind(sequence_info, batch_sequence_info)
        
        # Print progress
        cat("Retrieved sequences for batch:", start, "to", end, "on attempt", attempt, "\n")
        success <- TRUE
      }, silent = TRUE)
      
      if (!success) {
        cat("Failed to retrieve batch:", start, "to", end, "on attempt", attempt, "\n")
        attempt <- attempt + 1
      }
    }
    
    if (!success) {
      cat("Failed to retrieve batch:", start, "to", end, "after", max_retries, "attempts\n")
    }
  }
  
  return(sequence_info)
}

# Retrieve sequence information in batches
sequence_info <- get_sequences_in_batches(target_ensembl_ids, batch_size = 50)

# Display the retrieved sequence information
cat("Sequence information:\n")
print(head(sequence_info))

# Ensure target_ensembl in mirna_targets is character
mirna_targets$target_ensembl <- as.character(mirna_targets$target_ensembl)

# Ensure ensembl_gene_id in sequence_info is character
sequence_info$ensembl_gene_id <- as.character(sequence_info$ensembl_gene_id)

# Print the data types to ensure they match
cat("Structure of mirna_targets:\n")
str(mirna_targets)
cat("Structure of sequence_info:\n")
str(sequence_info)

# Print first few rows of mirna_targets and sequence_info before merging
cat("First few rows of mirna_targets:\n")
print(head(mirna_targets))

cat("First few rows of sequence_info:\n")
print(head(sequence_info))

# Print unique target_ensembl values in both data frames
cat("Unique target_ensembl in mirna_targets:\n")
print(unique(mirna_targets$target_ensembl))

cat("Unique ensembl_gene_id in sequence_info:\n")
print(unique(sequence_info$ensembl_gene_id))

# Merge the sequence information with mirna_targets
mirna_targets <- merge(mirna_targets, sequence_info, by.x = "target_ensembl", by.y = "ensembl_gene_id", all.x = TRUE)

# Display the merged data frame
cat("Merged data frame with sequences:\n")
print(head(mirna_targets))

# Write the merged data to a new CSV file
write.csv(mirna_targets, "/Volumes/WDT_GDS/mirna_targets_validated.csv", row.names = FALSE)