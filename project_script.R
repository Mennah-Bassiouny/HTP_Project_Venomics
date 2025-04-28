# Load required package
library(mzID)

# Step 1: Read the .mzid file
mzid_file <- "peptides_1_1_0_Ar_proteomics.mzid.xml"
mzid_data <- mzID(mzid_file)

# Step 2: Flatten the mzid object
mzid_flat <- flatten(mzid_data)

# Step 3: View column names to confirm structure
cat("🧬 Column names detected:\n")
print(colnames(mzid_flat))

# Step 4: Extract likely peptide + protein info columns
# This will include any columns with "pep", "acc", "charge", "mass", etc.
peptide_table <- mzid_flat[, grep("pep|acc|charge|mass|spectrum", names(mzid_flat), ignore.case = TRUE)]

# Step 5: View a sample of the table
head(peptide_table)

# Step 6: Save to CSV for downstream use
write.csv(peptide_table, "mzid_peptide_matches.csv", row.names = FALSE)
cat("✅ Exported: mzid_peptide_matches.csv\n")

################################################################################

# Load libraries
library(mzID)
library(readr)
library(dplyr)

# Step 1: Load mzID and flatten it
mzid_file <- "peptides_1_1_0_Ar_proteomics.mzid.xml"
mzid_data <- mzID(mzid_file)
mzid_flat <- flatten(mzid_data)

# Step 2: Load curated Table 2 proteins
table2 <- read_csv("venom_table2_curated_proteins.csv")

# Step 3: Check actual column names in mzid_flat
colnames(mzid_flat)

# Step 4: Identify the correct column with accessions (adjust if needed)
# Common ones: "accession", "accession.number", etc.
# You may need to update the column name here based on previous step
matched_peptides <- mzid_flat %>%
  filter(accession %in% table2$Accession)

# Step 5: Merge details from Table 2 (e.g., description, %)
matched_full <- matched_peptides %>%
  left_join(table2, by = c("accession" = "Accession"))

# Step 6: Save the result
write.csv(matched_full, "matched_peptides_to_table2.csv", row.names = FALSE)
cat("✅ Done: matched peptides written to matched_peptides_to_table2.csv\n")

################################################################################
