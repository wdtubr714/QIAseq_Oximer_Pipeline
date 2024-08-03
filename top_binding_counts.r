library(Biostrings)
library(dplyr)

# Set base directory
base_dir <- "/data4/msc19104442/target_ID"

# Load data
mirna_DE <- read.csv(file.path(base_dir, "mirna_DE.csv"))
mirna_targets_validated <- read.csv(file.path(base_dir, "mirna_targets_validated.csv"))
mirna_targets_predicted <- read.csv(file.path(base_dir, "mirna_targets_predicted.csv"))

# Print column names to ensure correct columns are being used
cat("Column names in mirna_targets_validated:\n")
print(colnames(mirna_targets_validated))

cat("Column names in mirna_DE:\n")
print(colnames(mirna_DE))

cat("Column names in mirna_targets_predicted:\n")
print(colnames(mirna_targets_predicted))

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
  
  # Check for the presence of 'U' or 'T' at the first position for "7mer-A1" and "8mer"
  if ((seed_type == "7mer_A1" || seed_type == "8mer") && !(startsWith(as.character(seed_seq), "U") || startsWith(as.character(seed_seq), "T"))) {
    return(0)
  }
  
  # Count matches of seed sequence reverse complement in the 3' UTR sequence
  binding_counts <- countPattern(seed_seq_rc, utr_seq, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
  
  return(binding_counts)
}

# Function to process data and generate top binding probability files
process_data <- function(targets, de_data, output_file_normal, output_file_mutated) {
  # Merge datasets using a left join to retain all target interactions
  merged_df <- merge(targets, de_data, by.x = "mature_mirna_id", by.y = "miRNA", all.x = TRUE)
  print(head(merged_df))
  
  # Filter out rows with NA values or 'Sequence unavailable' in the X3utr column
  filtered_df <- merged_df %>%
    filter(!is.na(X3utr) & X3utr != "Sequence unavailable")
  print(head(filtered_df))
  
  # Calculate binding scores for normal and mutated seed sequences
  seed_types <- c("6mer", "7mer_m8", "7mer_A1", "8mer")
  
  for (seed_type in seed_types) {
    filtered_df[[paste0("binding_score_normal_", seed_type)]] <- mapply(calculate_binding_score, 
                                                                        filtered_df$seed_sequence_normal, 
                                                                        filtered_df$X3utr, 
                                                                        MoreArgs = list(seed_type = seed_type))
    
    filtered_df[[paste0("binding_score_mutated_", seed_type)]] <- mapply(calculate_binding_score, 
                                                                         filtered_df$seed_sequence_mutated, 
                                                                         filtered_df$X3utr, 
                                                                         MoreArgs = list(seed_type = seed_type))
  }
  
  # Calculate total binding counts for normal and mutated sequences
  filtered_df <- filtered_df %>%
    mutate(total_binding_counts_normal = binding_score_normal_6mer + binding_score_normal_7mer_m8 + binding_score_normal_7mer_A1 + binding_score_normal_8mer,
           total_binding_counts_mutated = binding_score_mutated_6mer + binding_score_mutated_7mer_m8 + binding_score_mutated_7mer_A1 + binding_score_mutated_8mer)
  
  # Filter to keep rows with non-zero binding scores for normal sequences
  filtered_binding_data_normal <- filtered_df %>%
    filter(total_binding_counts_normal > 0)
  
  # Filter to keep rows with non-zero binding scores for mutated sequences
  filtered_binding_data_mutated <- filtered_df %>%
    filter(total_binding_counts_mutated > 0)
  
  # Get the top five binding probabilities for normal sequences
  top_five_binding_normal <- filtered_binding_data_normal %>%
    group_by(seed_sequence_normal) %>%
    arrange(desc(total_binding_counts_normal)) %>%  # Arrange in descending order
    slice_head(n = 5) %>%  # Select top 5 rows from each group
    ungroup() %>%  # Ungroup the data
    select(mature_mirna_id, pos.mut, database, target_ensembl, target_symbol, seed_sequence_normal, seed_sequence_mutated,
           binding_score_normal_6mer, binding_score_normal_7mer_m8, binding_score_normal_7mer_A1, binding_score_normal_8mer, total_binding_counts_normal)  # Select relevant columns
  
  # Get the top five binding probabilities for mutated sequences
  top_five_binding_mutated <- filtered_binding_data_mutated %>%
    group_by(seed_sequence_mutated) %>%
    arrange(desc(total_binding_counts_mutated)) %>%  # Arrange in descending order
    slice_head(n = 5) %>%  # Select top 5 rows from each group
    ungroup() %>%  # Ungroup the data
    select(mature_mirna_id, pos.mut, database, target_ensembl, target_symbol, seed_sequence_normal, seed_sequence_mutated,
           binding_score_mutated_6mer, binding_score_mutated_7mer_m8, binding_score_mutated_7mer_A1, binding_score_mutated_8mer, total_binding_counts_mutated)  # Select relevant columns
  
  # Write the results to a new CSV file
  write.csv(top_five_binding_normal, output_file_normal, row.names = FALSE)
  cat("Top five binding probabilities for normal sequences written to:", output_file_normal, "\n")
  
  write.csv(top_five_binding_mutated, output_file_mutated, row.names = FALSE)
  cat("Top five binding probabilities for mutated sequences written to:", output_file_mutated, "\n")
}

# Process validated targets
output_file_validated_normal <- file.path(base_dir, "top_five_binding_validated_normal.csv")
output_file_validated_mutated <- file.path(base_dir, "top_five_binding_validated_mutated.csv")
process_data(mirna_targets_validated, mirna_DE, output_file_validated_normal, output_file_validated_mutated)

# Process predicted targets
output_file_predicted_normal <- file.path(base_dir, "top_five_binding_predicted_normal.csv")
output_file_predicted_mutated <- file.path(base_dir, "top_five_binding_predicted_mutated.csv")
process_data(mirna_targets_predicted, mirna_DE, output_file_predicted_normal, output_file_predicted_mutated)