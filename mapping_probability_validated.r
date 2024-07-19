library(Biostrings)
library(dplyr)

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

# Extract miRNA ID from miRNA_MIMAT_ID_Species for merging
mirna_DE <- mirna_DE %>%
  mutate(miRNA_ID = sapply(strsplit(as.character(miRNA_MIMAT_ID_Species), " "), `[`, 1))
print(head(mirna_DE))

# Merge datasets using a left join to retain all target interactions
merged_df <- merge(mirna_targets, mirna_DE, by.x = "mature_mirna_id", by.y = "miRNA_ID", all.x = TRUE)
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

# Calculate binding scores for each seed type
for (type in seed_types) {
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(
      !!paste0("binding_score_normal_", gsub("-", "_", type)) := calculate_binding_score(seed_sequence_normal, X3utr, type),
      !!paste0("binding_score_mutated_", gsub("-", "_", type)) := calculate_binding_score(seed_sequence_mutated, X3utr, type)
    ) %>%
    ungroup()
}

# Function to normalize scores using Z-score within each miRNA group
z_score_normalization_group <- function(df, score_col) {
  df %>%
    group_by(mature_mirna_id) %>%
    mutate(
      mean_score = mean(!!sym(score_col), na.rm = TRUE),
      sd_score = sd(!!sym(score_col), na.rm = TRUE),
      z_score = ifelse(sd_score == 0, 0, (!!sym(score_col) - mean_score) / sd_score)
    ) %>%
    ungroup() %>%
    select(-mean_score, -sd_score)
}

# Normalize scores using Z-score for each seed type within miRNA groups
for (type in seed_types) {
  normal_col <- paste0("binding_score_normal_", gsub("-", "_", type))
  mutated_col <- paste0("binding_score_mutated_", gsub("-", "_", type))
  
  merged_df <- z_score_normalization_group(merged_df, normal_col) %>%
    rename(!!paste0("binding_prob_normal_", gsub("-", "_", type)) := z_score)
  
  merged_df <- z_score_normalization_group(merged_df, mutated_col) %>%
    rename(!!paste0("binding_prob_mutated_", gsub("-", "_", type)) := z_score)
}

# Summarize binding scores into a single combined score
merged_df <- merged_df %>%
  rowwise() %>%
  mutate(
    combined_binding_prob_normal = sum(abs(c_across(starts_with("binding_prob_normal_"))), na.rm = TRUE),
    combined_binding_prob_mutated = sum(abs(c_across(starts_with("binding_prob_mutated_"))), na.rm = TRUE)
  ) %>%
  ungroup()

# Add columns that count all the counts up for normal and mutated binding scores
merged_df <- merged_df %>%
  rowwise() %>%
  mutate(
    total_binding_counts_normal = sum(c_across(starts_with("binding_score_normal_")), na.rm = TRUE),
    total_binding_counts_mutated = sum(c_across(starts_with("binding_score_mutated_")), na.rm = TRUE)
  ) %>%
  ungroup()

# Print column names to verify the presence of all required columns
cat("Column names in merged_df:\n")
print(colnames(merged_df))

cat("miRNA determination completed.\n")

# Select final data frame including 3' UTR and 5' seed sequences and probabilities from Biostrings
cat("Selecting final data frame...\n")
final_df <- merged_df %>%
  select(miRNA_MIMAT_ID_Species, database, target_ensembl, pubmed_id, target_symbol, X3utr, adjusted_strand_with_mutation,
         seed_sequence_normal, seed_sequence_mutated, pos.mut, experiment,
         binding_prob_normal_8mer, binding_prob_normal_7mer_m8, binding_prob_normal_7mer_A1, binding_prob_normal_6mer,
         binding_prob_mutated_8mer, binding_prob_mutated_7mer_m8, binding_prob_mutated_7mer_A1, binding_prob_mutated_6mer,
         combined_binding_prob_normal, combined_binding_prob_mutated,
         total_binding_counts_normal, total_binding_counts_mutated,
         binding_score_normal_6mer, binding_score_normal_7mer_m8, binding_score_normal_7mer_A1, binding_score_normal_8mer,
         binding_score_mutated_6mer, binding_score_mutated_7mer_m8, binding_score_mutated_7mer_A1, binding_score_mutated_8mer)

cat("Final data frame selected.\n")

# Define output file path
output_file <- file.path(base_dir, "mirna_binding_probabilities_biostrings_validated.csv")

# Write the final data frame to a CSV file in the base directory
cat("Writing final data frame to CSV...\n")
write.csv(final_df, output_file, row.names = FALSE)

# Print message indicating that the file has been written
cat("Output file created at:", output_file, "\n")
