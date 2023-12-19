# Load the biomaRt package
library(biomaRt)
library(readxl)




n
# Define the Ensembl dataset and organism for C. elegans
ensembl_dataset <- "celegans_gene_ensembl"
ensembl_organism <- "celegans"

# Load the Ensembl database

mart <- useMart(biomart = "ensembl", dataset = ensembl_dataset)


# Define a list of common gene names (replace with your list)

OverlapAnal <- read_excel("C:/Users/jacob/School/Lab/CLIPvMS/OverlapAnal.xlsx")

wb_gene_ids <- getBM(
  attributes = c("ensembl_gene_id"),
  filters = "external_gene_name",
  values = OverlapAnal$GeneName,
  mart = mart
)

# Add the WormBase Gene IDs as a second column in your dataframe
OverlapAnal$GeneID <- wb_gene_ids$ensembl_gene_id

# Load the Ensembl database
ensembl <- useEnsembl(ensembl = ensembl_dataset, mart = ensembl_organism)

# Fetch the amino acid sequences
amino_acid_sequences <- getBM(
  attributes = c("ensembl_gene_id", "peptide"),
  filters = "wormbase_gene",
  values = OverlapAnal$GeneID,
  mart = mart
)

# Merge the amino acid sequences with your dataframe

merged_dataframe <- merge(amino_acid_sequences, OverlapAnal, by.x = "ensembl_gene_id", by.y = "GeneID", all.x = TRUE)
merged_dataframe$peptide <- substr(merged_dataframe$peptide, 1, nchar(merged_dataframe$peptide) - 1)

# Define the output file name
fasta_file <- "output.fasta"

# Open the file for writing
file_conn <- file(fasta_file, "w")

# Write each row as a FASTA sequence
for (i in 1:nrow(merged_dataframe)) {
  gene_name <- merged_dataframe$GeneName[i]
  sequence <- merged_dataframe$peptide[i]
  
  # Write the description and sequence to the file
  writeLines(paste(">", gene_name, sep = ""), file_conn)
  writeLines(sequence, file_conn)
}

# Close the file connection
close(file_conn)
cat("FASTA file", fasta_file, "created.\n")
