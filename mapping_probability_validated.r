library(Biostrings)
library(dplyr)
library(purrr)

# Set base directory
base_dir <- "/path/to/dir"

# Load data
mirna_DE <- read.csv(file.path(base_dir, "mirna_DE.csv"))
mirna_targets <- read.csv(file.path(base_dir, "mirna_targets_validated.csv"))

# Print column names to ensure correct columns are being used
cat("Column names in mirna_targets:\n")
print(colnames(mirna_targets))

cat("Column names in mirna_DE:\n")
print(colnames(mirna_DE))

# Merge datasets using a left join to retain all target interactions
merged_df <- merge(mirna_targets, mirna_DE, by.x = "mature_mirna_id", by.y = "miRNA", all.x = TRUE)
print(head(merged_df))

# Filter out rows with NA values or 'Sequence unavailable' in the X3utr column
merged_df <- merged_df %>%
  filter(!is.na(X3utr) & X3utr != "Sequence unavailable")
print(head(merged_df))

# Function to calculate binding score for different seed lengths and types
calculate_binding_score <- function(seed_seq, utr_seq, seed_type) {
  if (is.na(seed_seq) | is.na(utr_seq)) {
    return(0)
  }
  
  # Truncate seed sequence based on the seed type
  truncate_seed <- function(seq, type) {
    switch(type,
           "6mer" = substr(seq, 2, 7),
           "7mer_m8" = substr(seq, 2, 8),
           "7mer_A1" = substr(seq, 1, 7),
           "8mer" = substr(seq, 1, 8),
           seq)  # Default case returns the original sequence
  }
  
  seed_seq <- truncate_seed(seed_seq, seed_type)
  
  # Ensure sequences are in DNA format
  seed_seq_dna <- DNAString(seed_seq)
  seed_seq_rc <- as.character(reverseComplement(seed_seq_dna))
  
  # Check for the presence of an 'A' at the appropriate position in the reverse complement for 7mer-A1 and 8mer
  if ((seed_type == "7mer_A1" || seed_type == "8mer") && substr(seed_seq_rc, 1, 1) != "A") {
    return(0)
  }
  
  # Count matches of seed sequence reverse complement in the 3' UTR sequence
  binding_counts <- vcountPattern(seed_seq_rc, utr_seq, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
  
  return(binding_counts)
}

# Define seed types
seed_types <- c("6mer", "7mer_m8", "7mer_A1", "8mer")

# Function to calculate binding scores and probabilities for each seed type for each row
calculate_row_probabilities <- function(row, seed_types) {
  utr_length <- nchar(row$X3utr)
  
  normal_scores <- sapply(seed_types, function(type) {
    score <- calculate_binding_score(row$seed_sequence_normal, row$X3utr, type)
    return(score / utr_length)  # Normalize by UTR length
  })
  
  mutated_scores <- sapply(seed_types, function(type) {
    score <- calculate_binding_score(row$seed_sequence_mutated, row$X3utr, type)
    return(score / utr_length)  # Normalize by UTR length
  })
  
  # Sum normalized scores to get collective binding scores
  collective_normal_score <- sum(normal_scores)
  collective_mutated_score <- sum(mutated_scores)
  
  # Calculate probabilities
  probabilities <- sapply(1:length(seed_types), function(i) {
    normal_score <- normal_scores[i]
    mutated_score <- mutated_scores[i]
    if (normal_score == 0) {
      return(0)
    } else {
      return(1 - (mutated_score / normal_score))
    }
  })
  
  names(probabilities) <- paste0("mutation_prob_", seed_types)
  
  # Calculate binding counts for each seed type without normalization
  normal_binding_counts <- sapply(seed_types, function(type) {
    calculate_binding_score(row$seed_sequence_normal, row$X3utr, type)
  })
  
  mutated_binding_counts <- sapply(seed_types, function(type) {
    calculate_binding_score(row$seed_sequence_mutated, row$X3utr, type)
  })
  
  # Sum binding counts to get total counts
  total_normal_binding_counts <- sum(normal_binding_counts)
  total_mutated_binding_counts <- sum(mutated_binding_counts)
  
  binding_counts <- c(normal_binding_counts, mutated_binding_counts)
  names(binding_counts) <- c(paste0("normal_binding_count_", seed_types), paste0("mutated_binding_count_", seed_types))
  
  return(c(probabilities, collective_normal_score = collective_normal_score, collective_mutated_score = collective_mutated_score,
           total_normal_binding_counts = total_normal_binding_counts, total_mutated_binding_counts = total_mutated_binding_counts,
           binding_counts))
}

# Apply the function to each row in the dataframe using pmap
probabilities_list <- pmap(merged_df, calculate_row_probabilities, seed_types = seed_types)
probabilities_df <- as.data.frame(do.call(rbind, probabilities_list))
merged_df <- cbind(merged_df, probabilities_df)

# Function to extract top 3 target genes based on collective binding scores
extract_top_targets <- function(df, score_col, top_n = 3) {
  df %>%
    group_by(mature_mirna_id, pos.mut) %>%
    top_n(top_n, !!sym(score_col)) %>%
    ungroup()
}

# Extract top 3 targets for normal and mutated sequences
top_targets_normal <- extract_top_targets(merged_df, "collective_normal_score")
top_targets_mutated <- extract_top_targets(merged_df, "collective_mutated_score")

# Combine top targets into one dataframe
top_targets_combined <- bind_rows(top_targets_normal, top_targets_mutated)

# Print column names to verify the presence of all required columns
cat("Column names in top_targets_combined:\n")
print(colnames(top_targets_combined))

cat("miRNA determination completed.\n")

# Select final data frame including 3' UTR and 5' seed sequences and probabilities from Biostrings
cat("Selecting final data frame...\n")
final_df <- top_targets_combined %>%
  select(mature_mirna_id, database, target_ensembl, pubmed_id, target_symbol, X3utr, Sequence,
         seed_sequence_normal, seed_sequence_mutated, pos.mut, experiment,
         starts_with("mutation_prob_"), collective_normal_score, collective_mutated_score,
         total_normal_binding_counts, total_mutated_binding_counts,
         starts_with("normal_binding_count_"), starts_with("mutated_binding_count_"))

cat("Final data frame selected.\n")

# Define output file path
output_file <- file.path(base_dir, "mirna_binding_probabilities_biostrings_validated_top_targets.csv")

# Write the final data frame to a CSV file in the base directory
cat("Writing final data frame to CSV...\n")
write.csv(final_df, output_file, row.names = FALSE)

# Print message indicating that the file has been written
cat("Output file created at:", output_file, "\n")
