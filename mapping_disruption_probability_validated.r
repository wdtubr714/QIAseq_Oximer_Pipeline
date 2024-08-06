# Load required libraries
library(Biostrings)
library(dplyr)
library(purrr)

# Set base directory
base_dir <- "/path/to/dir"

# Load data
mirna_DE <- read.csv(file.path(base_dir, "mirna_DE.csv"))
mirna_targets_validated <- read.csv(file.path(base_dir, "mirna_targets_validated.csv"))

# Print column names to ensure correct columns are being used
cat("Column names in mirna_targets_validated:\n")
print(colnames(mirna_targets_validated))

cat("Column names in mirna_DE:\n")
print(colnames(mirna_DE))

# Merge datasets using a left join to retain all target interactions
merged_df <- merge(mirna_targets_validated, mirna_DE, by.x = "mature_mirna_id", by.y = "miRNA", all.x = TRUE)
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
  binding_counts <- countPattern(seed_seq_rc, utr_seq, max.mismatch=0, min.mismatch=0, with.indels=FALSE, fixed=TRUE)
  
  return(binding_counts)
}

# Define seed types
seed_types <- c("6mer", "7mer_m8", "7mer_A1", "8mer")

# Function to count possible binding sites in the UTR
count_possible_binding_sites <- function(utr_seq, seed_type) {
  # Define seed lengths
  seed_lengths <- c("6mer" = 6, "7mer_m8" = 7, "7mer_A1" = 7, "8mer" = 8)
  seed_length <- seed_lengths[seed_type]
  
  # Count possible binding sites
  possible_sites <- nchar(utr_seq) - seed_length + 1
  utr_length <- nchar(utr_seq)
  
  # Calculate density of binding sites
  binding_site_density <- possible_sites / utr_length
  
  return(list(possible_sites = possible_sites, binding_site_density = binding_site_density))
}

# Function to calculate binding scores and disruption probabilities for each seed type for each row
calculate_row_probabilities <- function(row, seed_types) {
  normal_scores <- sapply(seed_types, function(type) {
    calculate_binding_score(row$seed_sequence_normal, row$X3utr, type)
  })
  
  mutated_scores <- sapply(seed_types, function(type) {
    calculate_binding_score(row$seed_sequence_mutated, row$X3utr, type)
  })
  
  possible_sites_info <- lapply(seed_types, function(type) {
    count_possible_binding_sites(row$X3utr, type)
  })
  
  possible_sites <- sapply(possible_sites_info, function(info) info$possible_sites)
  binding_site_density <- sapply(possible_sites_info, function(info) info$binding_site_density)
  
  # Handle NA values in possible_sites and binding_site_density
  possible_sites[is.na(possible_sites)] <- 1  # To avoid division by zero
  binding_site_density[is.na(binding_site_density)] <- 1
  
  # Combine absolute counts and density for normalization
  normalization_factor <- (possible_sites + binding_site_density) / 2
  
  # Normalize scores based on the combined normalization factor
  normalized_normal_scores <- normal_scores / normalization_factor
  normalized_normal_scores[is.na(normalized_normal_scores)] <- 0  # Handle NA values
  
  normalized_mutated_scores <- mutated_scores / normalization_factor
  normalized_mutated_scores[is.na(normalized_mutated_scores)] <- 0  # Handle NA values
  
  # Calculate disruption probabilities for each seed type
  disruption_probabilities <- sapply(1:length(seed_types), function(i) {
    normalized_mutated_score <- normalized_mutated_scores[i]
    normalized_normal_score <- normalized_normal_scores[i]
    prob <- normalized_mutated_score / (normalized_mutated_score + normalized_normal_score + 1e-6)
    return(prob)
  })
  
  names(disruption_probabilities) <- paste0("disruption_prob_", seed_types)
  
  binding_counts <- c(normal_scores, mutated_scores)
  names(binding_counts) <- c(paste0("normal_binding_count_", seed_types), paste0("mutated_binding_count_", seed_types))
  
  # Calculate the collective disruption probability (average disruption probability)
  collective_disruption_probability <- mean(disruption_probabilities, na.rm = TRUE)
  
  return(c(disruption_probabilities, collective_normal_score = sum(normalized_normal_scores, na.rm = TRUE), 
           collective_mutated_score = sum(normalized_mutated_scores, na.rm = TRUE),
           collective_disruption_probability = collective_disruption_probability, binding_counts,
           total_normal_binding_counts = sum(normal_scores, na.rm = TRUE), 
           total_mutated_binding_counts = sum(mutated_scores, na.rm = TRUE)))
}

# Apply the function to each row in the dataframe using pmap
probabilities_list <- pmap(merged_df, ~ calculate_row_probabilities(list(...), seed_types = seed_types))
probabilities_df <- as.data.frame(do.call(rbind, probabilities_list))
merged_df <- cbind(merged_df, probabilities_df)

# Function to extract top 3 target genes based on collective disruption probability
extract_top_targets <- function(df, score_col, top_n = 3) {
  df %>%
    group_by(mature_mirna_id, pos.mut) %>%
    top_n(top_n, !!sym(score_col)) %>%
    ungroup()
}

# Extract top 3 targets for normal and mutated sequences
top_targets <- extract_top_targets(merged_df, "collective_disruption_probability")

# Print column names to verify the presence of all required columns
cat("Column names in top_targets:\n")
print(colnames(top_targets))

cat("miRNA determination completed.\n")

# Select final data frame including 3' UTR and 5' seed sequences and probabilities from Biostrings
cat("Selecting final data frame...\n")
final_df <- top_targets %>%
  select(mature_mirna_id, database, target_ensembl, pubmed_id, experiment, target_symbol, Sequence,
         seed_sequence_normal, seed_sequence_mutated, pos.mut,
         starts_with("disruption_prob_"), collective_normal_score, collective_mutated_score,
         collective_disruption_probability, starts_with("normal_binding_count_"), starts_with("mutated_binding_count_"),
         total_normal_binding_counts, total_mutated_binding_counts)

cat("Final data frame selected.\n")

# Define output file path
output_file <- file.path(base_dir, "mirna_binding_probabilities_disruption_validated_top_targets.csv")

# Write the final data frame to a CSV file in the base directory
cat("Writing final data frame to CSV...\n")
write.csv(final_df, output_file, row.names = FALSE)

# Print message indicating that the file has been written
cat("Output file created at:", output_file, "\n")
