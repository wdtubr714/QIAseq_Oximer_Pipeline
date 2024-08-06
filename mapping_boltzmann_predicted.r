# Load necessary libraries
library(Biostrings)
library(dplyr)
library(purrr)
library(pbapply)

# Set base directory
base_dir <- "/path/to/dir"

# Load data
mirna_DE <- read.csv(file.path(base_dir, "mirna_DE_oxoG.csv"))
mirna_targets_predicted <- read.csv(file.path(base_dir, "mirna_targets_predicted.csv"))

# Print column names to ensure correct columns are being used
cat("Column names in mirna_targets_predicted:\n")
print(colnames(mirna_targets_predicted))

cat("Column names in mirna_DE:\n")
print(colnames(mirna_DE))

# Merge datasets using a left join to retain all target interactions
merged_df <- merge(mirna_targets_predicted, mirna_DE, by.x = "mature_mirna_id", by.y = "miRNA", all.x = TRUE)
print(head(merged_df))

# Filter out rows with NA values or 'Sequence unavailable' in the X3utr column
merged_df <- merged_df %>%
  filter(!is.na(X3utr) & X3utr != "Sequence unavailable")
print(head(merged_df))

# Filter out rows with NA values or 'Sequence unavailable' in the pos.mut column
merged_df <- merged_df %>%
  filter(!is.na(pos.mut))
print(head(merged_df))

# Function to truncate seed sequence based on its type
truncate_seed <- function(seq, type) {
  switch(type,
         "6mer" = substr(seq, 2, 7),
         "7mer_m8" = substr(seq, 2, 8),
         "7mer_A1" = substr(seq, 1, 7),
         "8mer" = substr(seq, 1, 8),
         seq)  # Default case returns the original sequence
}

calculate_binding_score <- function(seed_seq, utr_seq, seed_type) {
  if (is.na(seed_seq) | is.na(utr_seq)) {
    return(0)
  }
  
  # Truncate seed sequence based on the seed type
  seed_seq <- truncate_seed(seed_seq, seed_type)
  
  # Convert 'T' to 'U' for RNA sequences
  seed_seq_rna <- chartr("T", "U", seed_seq)
  
  # Ensure sequences are in RNA format
  seed_seq_rna <- RNAString(seed_seq_rna)
  seed_seq_rc <- as.character(reverseComplement(seed_seq_rna))
  
  # Convert utr_seq to RNA
  utr_seq_rna <- chartr("T", "U", utr_seq)
  utr_seq_rna <- RNAString(utr_seq_rna)
  
  # Check for the presence of 'U' or 'T' at the first position for "7mer-A1" and "8mer"
  if ((seed_type == "7mer_A1" || seed_type == "8mer") && !(startsWith(as.character(seed_seq_rna), "U") || startsWith(as.character(seed_seq_rna), "T"))) {
    cat("Caution: Seed sequence type", seed_type, "does not meet the criteria and will have zero binding count.\n")
    return(0)
  }
  
  # Count matches of seed sequence reverse complement in the 3' UTR sequence
  binding_counts <- countPattern(seed_seq_rc, utr_seq_rna, max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, fixed = TRUE)
  
  return(binding_counts)
}

# Function to calculate Gibbs free energy for a given seed type
calculate_gibbs_free_energy_for_seed <- function(seed_seq, utr_seq, seed_type) {
  # Truncate seed sequence based on seed type
  truncated_seed_seq <- truncate_seed(seed_seq, seed_type)
  
  # Convert sequences to RNA
  truncated_seed_seq <- chartr("T", "U", truncated_seed_seq)
  utr_seq <- chartr("T", "U", utr_seq)
  
  # Convert to RNAString
  truncated_seed_seq <- RNAString(truncated_seed_seq)
  utr_seq_rna <- RNAString(utr_seq)
  
  # Write sequences to temporary files
  input_file <- tempfile(pattern = "input", fileext = ".fa")
  output_file <- tempfile(pattern = "duplex", fileext = ".out")
  
  writeLines(c(paste(">seed", truncated_seed_seq, sep="\n"), paste(">utr", utr_seq, sep="\n")), input_file)
  
  # Print contents of the temporary input file for debugging
  cat("Contents of input file for seed type", seed_type, ":\n")
  print(readLines(input_file))
  
  # Use RNAduplex from ViennaRNA package to calculate free energy
  system(paste("RNAduplex < ", input_file, " > ", output_file, sep = ""))
  
  # Read the output and extract ΔG value
  duplex_out <- readLines(output_file)
  cat("RNAduplex output for seed type", seed_type, ":\n")
  print(duplex_out)
  
  delta_g <- NA  # Initialize delta_g as NA
  if (length(duplex_out) > 2 && grepl("\\(.*\\)", duplex_out[3])) {
    delta_g <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", duplex_out[3]))
    if (is.na(delta_g)) {
      warning("Failed to extract ΔG value for seed type", seed_type, ".")
    }
  } else {
    warning("RNAduplex output for seed type", seed_type, "is empty or does not contain expected format.")
  }
  cat("Extracted ΔG value for seed type", seed_type, ":", delta_g, "\n")
  
  # Clean up temporary files
  file.remove(input_file, output_file)
  
  return(delta_g)
}

# Function to calculate binding probability based on Gibbs free energy
calculate_binding_probability <- function(delta_g, temperature=310.15) {
  R <- 1.987  # Universal gas constant in cal/(K*mol)
  
  if (is.na(delta_g)) {
    return(NA)
  }
  
  # Calculate probability using the Boltzmann distribution
  p_binding <- exp(-delta_g / (R * temperature)) / (1 + exp(-delta_g / (R * temperature)))
  cat("Calculated binding probability for ΔG =", delta_g, "is", p_binding, "\n")
  
  return(p_binding)
}

# Function to calculate binding probabilities for each seed type
calculate_binding_probability_for_seed <- function(seed_seq, utr_seq, seed_type, temperature=310.15) {
  # Truncate seed sequence based on the seed type
  truncated_seed_seq <- truncate_seed(seed_seq, seed_type)
  
  # Convert sequences to RNA
  truncated_seed_seq <- chartr("T", "U", truncated_seed_seq)
  utr_seq <- chartr("T", "U", utr_seq)
  
  # Convert to RNAString
  truncated_seed_seq <- RNAString(truncated_seed_seq)
  utr_seq_rna <- RNAString(utr_seq)
  
  # Check for the presence of 'U' or 'T' at the first position for "7mer-A1" and "8mer"
  if ((seed_type == "7mer_A1" || seed_type == "8mer") && !(startsWith(as.character(truncated_seed_seq), "U") || startsWith(as.character(truncated_seed_seq), "T"))) {
    cat("Caution: Seed sequence type", seed_type, "does not meet the criteria and will have 0 probability.\n")
    return(0)
  }
  
  # Calculate Gibbs free energy
  delta_g <- calculate_gibbs_free_energy_for_seed(seed_seq, utr_seq, seed_type)
  
  # Calculate binding probability
  p_binding <- calculate_binding_probability(delta_g, temperature)
  
  return(p_binding)
}

calculate_row_probabilities <- function(row, seed_types) {
  # Print start message for the current row
  cat("Processing row with mature_mirna_id", row$mature_mirna_id, "...\n")
  
  # Initialize vectors to store scores and probabilities
  normal_scores <- numeric(length(seed_types))
  mutated_scores <- numeric(length(seed_types))
  probabilities_normal <- numeric(length(seed_types))
  probabilities_mutated <- numeric(length(seed_types))
  
  for (type in seed_types) {
    index <- which(seed_types == type)
    
    # Calculate binding scores
    normal_scores[index] <- calculate_binding_score(row$seed_sequence_normal, row$X3utr, type)
    mutated_scores[index] <- calculate_binding_score(row$seed_sequence_mutated, row$X3utr, type)
    
    # Calculate binding probabilities
    probabilities_normal[index] <- calculate_binding_probability_for_seed(row$seed_sequence_normal, row$X3utr, type)
    probabilities_mutated[index] <- calculate_binding_probability_for_seed(row$seed_sequence_mutated, row$X3utr, type)
  }
  
  # Filter out seed types with zero probabilities
  valid_indices <- which(probabilities_normal > 0 & probabilities_mutated > 0)
  
  # Filter scores and probabilities for valid seed types
  valid_seed_types <- seed_types[valid_indices]
  valid_normal_scores <- normal_scores[valid_indices]
  valid_mutated_scores <- mutated_scores[valid_indices]
  valid_probabilities_normal <- probabilities_normal[valid_indices]
  valid_probabilities_mutated <- probabilities_mutated[valid_indices]
  
  # Average the probabilities across the valid seed types
  avg_prob_normal <- if (length(valid_probabilities_normal) > 0) mean(valid_probabilities_normal, na.rm = TRUE) else NA
  avg_prob_mutated <- if (length(valid_probabilities_mutated) > 0) mean(valid_probabilities_mutated, na.rm = TRUE) else NA
  
  # Calculate collective scores and counts
  collective_normal_score <- sum(valid_normal_scores, na.rm = TRUE)
  collective_mutated_score <- sum(valid_mutated_scores, na.rm = TRUE)
  
  total_normal_binding_counts <- sum(valid_normal_scores, na.rm = TRUE)
  total_mutated_binding_counts <- sum(valid_mutated_scores, na.rm = TRUE)
  
  # Print final results for debugging
  cat("Final average binding probability for normal sequence:", avg_prob_normal, "\n")
  cat("Final average binding probability for mutated sequence:", avg_prob_mutated, "\n")
  cat("Collective normal score:", collective_normal_score, "\n")
  cat("Collective mutated score:", collective_mutated_score, "\n")
  
  # Return binding counts and probabilities
  return(c(avg_prob_normal = avg_prob_normal,
           avg_prob_mutated = avg_prob_mutated,
           collective_normal_score = collective_normal_score,
           collective_mutated_score = collective_mutated_score,
           total_normal_binding_counts = total_normal_binding_counts,
           total_mutated_binding_counts = total_mutated_binding_counts,
           setNames(as.list(normal_scores), paste0("normal_binding_count_", seed_types)),
           setNames(as.list(mutated_scores), paste0("mutated_binding_count_", seed_types))))
}

# Define seed types
seed_types <- c("6mer", "7mer_m8", "7mer_A1", "8mer")

# Apply the function to each row in the dataframe using pblapply
cat("Starting probability calculations for each row...\n")
probabilities_list <- pblapply(1:nrow(merged_df), function(i) {
  calculate_row_probabilities(merged_df[i, ], seed_types = seed_types)
})
cat("Probability calculations completed.\n")

# Combine results into a data frame
probabilities_df <- as.data.frame(do.call(rbind, probabilities_list))
merged_df <- cbind(merged_df, probabilities_df)

# Ensure avg_prob_normal and avg_prob_mutated are numeric
merged_df$avg_prob_normal <- as.numeric(as.character(merged_df$avg_prob_normal))
merged_df$avg_prob_mutated <- as.numeric(as.character(merged_df$avg_prob_mutated))

# Function to extract top n target genes based on binding probabilities
extract_top_targets <- function(df, prob_col, top_n = 3) {
  df %>%
    group_by(mature_mirna_id, pos.mut) %>%
    arrange(desc(!!sym(prob_col))) %>%
    slice_head(n = top_n) %>%
    ungroup()
}

# Extract top 3 targets for normal and mutated sequences
top_targets_normal <- extract_top_targets(merged_df, "avg_prob_normal")
top_targets_mutated <- extract_top_targets(merged_df, "avg_prob_mutated")

# Combine top targets into one dataframe
top_targets_combined <- bind_rows(top_targets_normal, top_targets_mutated)

# Print column names to verify the presence of all required columns
cat("Column names in top_targets_combined:\n")
print(colnames(top_targets_combined))

cat("miRNA determination completed.\n")

# Select final data frame including 3' UTR and 5' seed sequences and probabilities
cat("Selecting final data frame...\n")
final_df <- top_targets_combined %>%
  select(mature_mirna_id, database, target_ensembl, target_symbol, X3utr, Sequence,
         seed_sequence_normal, seed_sequence_mutated, pos.mut,
         avg_prob_normal, avg_prob_mutated, collective_normal_score, collective_mutated_score,
         starts_with("normal_binding_count_"), starts_with("mutated_binding_count_"))

# Convert list columns to character
final_df <- final_df %>%
  mutate(across(where(is.list), ~ sapply(., toString)))

cat("Final data frame selected.\n")

# Define output file path
output_file <- file.path(base_dir, "mirna_binding_probabilities_predicted_top_targets.csv")

# Write the final data frame to a CSV file in the base directory
cat("Writing final data frame to CSV...\n")
write.csv(final_df, output_file, row.names = FALSE)

# Print message indicating that the file has been written
cat("Output file created at:", output_file, "\n")

