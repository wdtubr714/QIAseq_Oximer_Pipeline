# Load necessary libraries
library(dplyr)

# Set base directory and define input and output file paths
base_dir <- "/data4/msc19104442/target_ID"
input_file <- file.path(base_dir, "mirna_binding_probabilities_biostrings_validated.csv")
mirna_de_file <- file.path(base_dir, "mirna_DE.csv")
output_file_normal <- file.path(base_dir, "top_five_binding_prob_normal_validated.csv")
output_file_mutated <- file.path(base_dir, "top_five_binding_prob_mutated_validated.csv")

# Load data
binding_data <- read.csv(input_file)
mirna_DE <- read.csv(mirna_de_file)

# Extract miRNA IDs from mirna_DE
mirna_ids <- unique(mirna_DE$miRNA)

# Filter binding data to include only miRNA types in mirna_DE
filtered_binding_data <- binding_data %>%
  filter(mature_mirna_id %in% mirna_ids)

# Get the top five binding probabilities for normal sequences
top_five_binding_normal <- filtered_binding_data %>%
  group_by(seed_sequence_normal) %>%
  arrange(desc(combined_binding_prob_normal)) %>%  # Arrange in descending order
  slice_head(n = 5) %>%  # Select top 5 rows from each group
  ungroup() %>%  # Ungroup the data
  select(mature_mirna_id, database, target_ensembl, pubmed_id, target_symbol, experiment,
         seed_sequence_normal, seed_sequence_mutated, pos.mut,
         binding_prob_normal_8mer, binding_prob_normal_7mer_m8, binding_prob_normal_7mer_A1, binding_prob_normal_6mer,
         combined_binding_prob_normal, binding_score_normal_6mer, binding_score_normal_7mer_m8,
         binding_score_normal_7mer_A1, binding_score_normal_8mer, total_binding_counts_normal)  # Select relevant columns

# Write the results to a new CSV file
write.csv(top_five_binding_normal, output_file_normal, row.names = FALSE)

# Print a message indicating that the file has been written
cat("Top five binding probabilities for normal sequences written to:", output_file_normal, "\n")

# Get the top five binding probabilities for mutated sequences
top_five_binding_mutated <- filtered_binding_data %>%
  group_by(seed_sequence_mutated) %>%
  arrange(desc(combined_binding_prob_mutated)) %>%  # Arrange in descending order
  slice_head(n = 5) %>%  # Select top 5 rows from each group
  ungroup() %>%  # Ungroup the data
  select(mature_mirna_id, database, target_ensembl, pubmed_id, target_symbol, experiment,
         seed_sequence_normal, seed_sequence_mutated, pos.mut,
         binding_prob_mutated_8mer, binding_prob_mutated_7mer_m8, binding_prob_mutated_7mer_A1, binding_prob_mutated_6mer,
         combined_binding_prob_mutated, binding_score_mutated_6mer, binding_score_mutated_7mer_m8,
         binding_score_mutated_7mer_A1, binding_score_mutated_8mer, total_binding_counts_mutated)  # Select relevant columns

# Write the results to a new CSV file
write.csv(top_five_binding_mutated, output_file_mutated, row.names = FALSE)

# Print a message indicating that the file has been written
cat("Top five binding probabilities for mutated sequences written to:", output_file_mutated, "\n")

