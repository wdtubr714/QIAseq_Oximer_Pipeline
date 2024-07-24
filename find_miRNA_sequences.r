# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# Directory containing the files
directory_path <- "/path/to/dir"

# List all .bowtie.mapped.txt files
file_list <- list.files(path = directory_path, pattern = "*.bowtie_mapped.txt", full.names = TRUE)

# Load miRNA list
miRNA_list_path <- file.path(directory_path, "miRNA_list.csv")
miRNA_list <- read_csv(miRNA_list_path)

# Assuming miRNA_list.csv has columns "miRNA" and "pos.mut"
colnames(miRNA_list) <- c("miRNA", "pos.mut")

# Convert miRNA names from miR to mir for consistency
miRNA_list <- miRNA_list %>%
  mutate(miRNA = str_replace_all(tolower(miRNA), "miR-", "mir-"))

# Function to convert 1-indexed positions to 0-indexed
convert_to_0_index <- function(pos_mut) {
  if (is.na(pos_mut) || pos_mut == "") return(pos_mut)
  
  mutations <- str_split(pos_mut, ",")[[1]]
  adjusted_mutations <- lapply(mutations, function(mut) {
    pos <- as.integer(str_extract(mut, "^\\d+")) - 1
    base_change <- str_extract(mut, ":[A-Z]>[A-Z]$")
    return(paste0(pos, base_change))
  })
  
  return(paste(adjusted_mutations, collapse = ","))
}

# Adjust positions in pos.mut to 0-index and store both 0-index and 1-index
miRNA_list <- miRNA_list %>%
  mutate(pos.mut_0_index = sapply(pos.mut, convert_to_0_index))

# Initialize a list to store results
all_results <- list()

# Process each file individually
for (file_path in file_list) {
  # Read the file
  data <- read_tsv(file_path, col_names = FALSE)
  
  # Add appropriate column names
  colnames(data) <- c("miRNA_MIMAT_ID_Species", "Direction", "Strand_with_count", "Count", "Sequence", "Quality", "Unknown", "Mutation", "miRNA")
  
  # Convert miRNA to lowercase for case-insensitive matching
  data <- data %>%
    mutate(miRNA_id = str_replace_all(tolower(str_extract(miRNA_MIMAT_ID_Species, "^[^ ]+")), "miR-", "mir-"))
  
  # Filter by both 0-index and 1-index mutations from miRNA_list
  data_filtered <- data %>%
    filter(Mutation %in% miRNA_list$pos.mut | Mutation %in% miRNA_list$pos.mut_0_index)
  
  # Filter by miRNA
  data_filtered <- data_filtered %>%
    filter(miRNA_id %in% miRNA_list$miRNA)
  
  # Log the data being processed
  cat("Processing file:", file_path, "\n")
  cat("Filtered data sample:\n")
  print(head(data_filtered))
  
  # Initialize a list to store file-specific results
  file_results <- list()
  
  # Process each miRNA individually
  for (i in 1:nrow(miRNA_list)) {
    miRNA <- miRNA_list$miRNA[i]
    pos_mut <- miRNA_list$pos.mut[i]
    pos_mut_0_index <- miRNA_list$pos.mut_0_index[i]
    
    # Filter the data for the specific miRNA and both mutation types
    filtered_data <- data_filtered %>%
      filter(miRNA_id == miRNA & (Mutation == pos_mut | Mutation == pos_mut_0_index))
    
    if (nrow(filtered_data) > 0) {
      # Select relevant columns
      filtered_data <- filtered_data %>%
        select(miRNA_MIMAT_ID_Species, miRNA_id, Strand_with_count, Mutation, Sequence)
      
      # Extract the numeric part from miRNA_MIMAT_ID_Species for sorting
      filtered_data <- filtered_data %>%
        mutate(sort_key = as.numeric(str_extract(miRNA_MIMAT_ID_Species, "(?<=miR-)[0-9]+")))
      
      # Sort the data based on the numeric part
      sorted_data <- filtered_data %>%
        arrange(sort_key) %>%
        select(-sort_key)
      
      # Append the results to the file-specific list
      file_results[[paste0(miRNA, "_", pos_mut)]] <- sorted_data
    }
  }
  
  # Combine the file-specific results into a single data frame
  if (length(file_results) > 0) {
    file_combined_data <- bind_rows(file_results)
    
    # Append the file-specific combined data to the overall results list
    all_results[[file_path]] <- file_combined_data
  }
}

# Combine all the results into a single data frame
final_data <- bind_rows(all_results)

# Log the final data before deduplication
cat("Final data sample before deduplication:\n")
print(head(final_data))

# Save the filtered final data to a file if needed
write_csv(final_data, file.path(directory_path, "filtered_final_sorted_filtered_miRNA_data.csv"))

# Deduplicate by matching miRNAs and their mutations from miRNA_list.csv
deduplicated_final_data <- final_data %>%
  inner_join(miRNA_list, by = c("miRNA_id" = "miRNA", "Mutation" = "pos.mut"), relationship = "many-to-many") %>%
  bind_rows(
    final_data %>%
      inner_join(miRNA_list, by = c("miRNA_id" = "miRNA", "Mutation" = "pos.mut_0_index"), relationship = "many-to-many")
  ) %>%
  distinct(Strand_with_count, .keep_all = TRUE)

# Ensure all 126 entries from miRNA_list are included
missing_entries <- miRNA_list %>%
  filter(!(miRNA %in% deduplicated_final_data$miRNA_id))

# Print missing entries if any
cat("Missing entries:\n")
print(missing_entries)

# Filter the deduplicated final data to only include the X miRNAs from the miRNA_list
final_filtered_data <- deduplicated_final_data %>%
  filter(miRNA_id %in% miRNA_list$miRNA)

# Save the final filtered data
write_csv(final_filtered_data, file.path(directory_path, "final_filtered_miRNA_data.csv"))

# Save the missing entries to a file if needed
write_csv(missing_entries, file.path(directory_path, "missing_miRNA_entries.csv"))

# Combine the final filtered data with miRNA_list to ensure only X miRNAs are included
combined_miRNA_list <- miRNA_list %>%
  inner_join(final_filtered_data %>% select(miRNA_id, Sequence), by = c("miRNA" = "miRNA_id"), relationship = "many-to-many")

# Reorder columns as miRNA, pos.mut, pos.mut_0_index, Sequence
combined_miRNA_list <- combined_miRNA_list %>%
  select(miRNA, pos.mut, pos.mut_0_index, Sequence) %>%
  distinct()

# Convert miRNA names from mir back to miR for consistency
combined_miRNA_list <- combined_miRNA_list %>%
  mutate(miRNA = str_replace_all(tolower(miRNA), "mir-", "miR-"))

# Extract the seed sequence (1 to 8) from the Sequence
combined_miRNA_list <- combined_miRNA_list %>%
  mutate(seed_sequence_mutated = substr(Sequence, 1, 8))

# Function to reverse mutations in a sequence
reverse_mutations <- function(sequence, mutations) {
  seq_list <- str_split(sequence, "")[[1]]
  if (!is.na(mutations)) {
    for (mutation in str_split(mutations, ",")[[1]]) {
      pos <- as.numeric(str_extract(mutation, "^\\d+"))
      base_change <- str_split(str_split(mutation, ":")[[1]][2], ">")[[1]]
      if (pos > 0 && pos <= length(seq_list)) {
        seq_list[pos] <- base_change[1]
      }
    }
  }
  return(paste(seq_list, collapse = ""))
}

# Create seed_sequence_normal by reversing the mutations
combined_miRNA_list <- combined_miRNA_list %>%
  rowwise() %>%
  mutate(seed_sequence_normal = reverse_mutations(seed_sequence_mutated, pos.mut))

# Save the updated miRNA_list with seed sequences
write_csv(combined_miRNA_list, file.path(directory_path, "mirna_DE.csv"))

# Log the combined miRNA_list
cat("Updated miRNA list with seed sequences sample:\n")
print(head(combined_miRNA_list))
