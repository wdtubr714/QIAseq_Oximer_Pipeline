# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)

# Directory containing the files
directory_path <- "/data4/msc19104442/O8G-miseq/"

# List all .bowtie.mapped.txt files
file_list <- list.files(path = directory_path, pattern = "*.bowtie.mapped.txt", full.names = TRUE)

# Read and merge all files
merged_data <- file_list %>% 
  lapply(read_tsv, col_names = FALSE) %>% 
  bind_rows()

# Add appropriate column names
colnames(merged_data) <- c("miRNA_MIMAT_ID_Species", "Direction", "Strand_with_count", "Count", "Sequence", "Quality", "Unknown", "Mutation", "miRNA")

# Extract the miRNA ID
merged_data <- merged_data %>%
  mutate(miRNA_id = str_extract(miRNA_MIMAT_ID_Species, "^[^ ]+"))

# Load miRNA list
miRNA_list_path <- file.path(directory_path, "miRNA_list.csv")
miRNA_list <- read_csv(miRNA_list_path)

# Assuming miRNA_list.csv has columns "miRNA" and "pos.mut"
colnames(miRNA_list) <- c("miRNA", "pos.mut")

# Convert miRNA_id to lowercase for case-insensitive matching
merged_data <- merged_data %>%
  mutate(miRNA_id = tolower(miRNA_id))

miRNA_list <- miRNA_list %>%
  mutate(miRNA = tolower(miRNA))

# Function to adjust positions in pos.mut
adjust_pos_mut <- function(pos_mut) {
  if (is.na(pos_mut) || pos_mut == "") return(pos_mut)
  
  mutations <- str_split(pos_mut, ",")[[1]]
  adjusted_mutations <- lapply(mutations, function(mut) {
    pos <- as.integer(str_extract(mut, "^\\d+")) - 1
    base_change <- str_extract(mut, ":[A-Z]>[A-Z]$")
    return(paste0(pos, base_change))
  })
  
  return(paste(adjusted_mutations, collapse = ","))
}

# Adjust positions in pos.mut
miRNA_list <- miRNA_list %>%
  mutate(pos.mut = sapply(pos.mut, adjust_pos_mut))

# Filter the merged data based on the miRNA list
filtered_data <- merged_data %>%
  inner_join(miRNA_list, by = c("miRNA_id" = "miRNA"), relationship = "many-to-many")

# Function to apply mutations to the sequence
apply_mutations <- function(sequence, pos_mut) {
  if (is.na(pos_mut) || pos_mut == "") return(sequence)
  
  mutations <- str_split(pos_mut, ",")[[1]]
  seq_split <- str_split(sequence, "")[[1]]
  
  for (mut in mutations) {
    pos <- as.integer(str_extract(mut, "^\\d+"))
    base <- str_extract(mut, "[A-Z]$")
    seq_split[pos + 1] <- base  # Adjusted for 1-based indexing in R
  }
  
  return(paste(seq_split, collapse = ""))
}

# Apply the function to the filtered data to create a new column for the mutated sequence
filtered_data <- filtered_data %>%
  mutate(
    adjusted_strand_with_mutation = mapply(apply_mutations, Sequence, pos.mut)
  ) %>%
  select(miRNA_MIMAT_ID_Species, Strand_with_count, pos.mut, adjusted_strand_with_mutation)

# Remove duplicate adjusted_strand_with_mutation
filtered_data <- filtered_data %>%
  distinct(adjusted_strand_with_mutation, .keep_all = TRUE)

# Extract the numeric part from miRNA_MIMAT_ID_Species for sorting
filtered_data <- filtered_data %>%
  mutate(sort_key = as.numeric(str_extract(miRNA_MIMAT_ID_Species, "(?<=miR-)[0-9]+")))

# Sort the data based on the numeric part
sorted_data <- filtered_data %>%
  arrange(sort_key) %>%
  select(-sort_key)

# Print the first few rows of the sorted data to confirm
print("Sorted data:")
print(head(sorted_data, 20))  # Print more rows for confirmation

# Save the sorted data to a file if needed
write_csv(sorted_data, file.path(directory_path, "sorted_filtered_miRNA_data.csv"))


