#!/usr/bin/env Rscript
library(data.table)
library(readr)

# Define the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Debugging: Print the received arguments
message("Received arguments: ", paste(args, collapse = ", "))

if (length(args) == 0) {
  stop("No input file provided. Please specify the TSV file as a command line argument.")
}

# Path to the input TSV file
tsv_file <- as.character(args[1])

# Print debugging information
message("Checking if TSV file exists: ", tsv_file)

# Check if the TSV file exists before attempting to read it
if (file.exists(tsv_file)) {
  message("TSV file exists: ", tsv_file)
  
  # Check the permissions of the file
  file_info <- file.info(tsv_file)
  print(file_info)
  
  # Print the first few lines of the file for inspection
  message("Printing the first few lines of the TSV file:")
  tryCatch({
    print(readLines(tsv_file, n = 5))
  }, error = function(e) {
    message("Error reading lines from TSV file: ", e$message)
  })
  
  # Read the TSV file using read_delim for better handling of delimiters
  message("Attempting to read the TSV file with read_delim:")
  tryCatch({
    df <- read_delim(tsv_file, delim = "\t")
    print(head(df))
    
    # Processing the data as per the provided logic
    message("Processing the data")
    counts <- as.data.table(df[!duplicated(df[["UID"]]), c(3, 13:ncol(df))])
    mirna <- counts[, lapply(.SD, sum), by = miRNA]
    
    # Write the output to a new file
    output_file <- file.path(dirname(tsv_file), "mirna.tsv")
    write.table(mirna, output_file, quote = FALSE, sep = "\t", row.names = FALSE)
    
    message("Output written to: ", output_file)
  }, error = function(e) {
    message("Error processing TSV file: ", e$message)
  })
} else {
  message("TSV file does not exist: ", tsv_file)
}

message("R script execution complete")