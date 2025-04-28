library(readr)
library(ggplot2)

# Step 1: Read the header line (5th line)
header_line <- readLines("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/eggNOG/gioj_annotation.emapper.annotations")[5]
col_names <- strsplit(gsub("^#", "", header_line), "\t")[[1]]

# Step 2: Read the data, skipping the first 4 lines and applying the column names
annotations <- read_tsv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/eggNOG/gioj_annotation.emapper.annotations", skip = 5, col_names = col_names)

# View the result
head(annotations)
#---------------------------------------

# View how GO and KEGG data look
head(annotations$GOs)
head(annotations$KEGG_ko)

# Optional: Count number of entries with GO or KEGG annotations
table(!is.na(annotations$GOs))   # TRUE means has GO annotation
table(!is.na(annotations$KEGG_ko))
#---------------------------------------

# Parse and clean the GO Terms
# Clean GO terms
go_terms <- annotations %>%
  filter(!is.na(GOs)) %>%
  separate_rows(GOs, sep = ",") %>%
  count(GOs, sort = TRUE)

# View top GO terms
head(go_terms, 10)
#---------------------------------------

# Parse KEGG Orthologs
# Count KEGG orthologs
kegg_ko_summary <- annotations %>%
  filter(!is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  count(KEGG_ko, sort = TRUE)

head(kegg_ko_summary, 10)
#---------------------------------------

# Visualize Top Terms
# Top 10 GO terms
go_terms %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(GOs, n), y = n)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  coord_flip() +
  labs(title = "Top 10 GO Terms", x = "GO Term", y = "Count") +
  theme_minimal()

# Top 10 KEGG KOs
kegg_ko_summary %>%
  slice_max(n, n = 10) %>%
  ggplot(aes(x = reorder(KEGG_ko, n), y = n)) +
  geom_bar(stat = "identity", fill = "#009E73") +
  coord_flip() +
  labs(title = "Top 10 KEGG Orthologs", x = "KEGG KO", y = "Count") +
  theme_minimal()
#-------------------------------------

# Install KEGGREST if not already installed
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}

library(KEGGREST)

# Example KO list
ko_list <- c("K11433", "K00699", "K06515", "K10455", "K09228", "K06767", "K10577", "K08202", "K15001")

# Function to retrieve KEGG pathway info for each KO
get_kegg_pathways <- function(ko_id) {
  tryCatch({
    pathway_links <- keggLink("pathway", ko_id)
    if (length(pathway_links) == 0) return(NULL)
    
    sapply(pathway_links, function(path_id) {
      keggGet(path_id)[[1]]$NAME
    })
  }, error = function(e) NA)
}

# Create dataframe
kegg_df <- data.frame(
  KO = ko_list,
  Pathways = sapply(ko_list, function(ko) {
    paste(get_kegg_pathways(ko), collapse = "; ")
  })
)

print(kegg_df)

# Save KEGG mappings to CSV
write.csv(kegg_df, "kegg_ko_pathway_mapping.csv", row.names = FALSE)

#------------------------------------

# Install required packages
if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")

library(httr)
library(jsonlite)

# Example GO terms
go_ids <- c("GO:0008150", "GO:0005575", "GO:0044464", "GO:0005623", "GO:0009987", "GO:0005622", "GO:0044424", "GO:0003674", "GO:0043226", "GO:0043229")

# Safe null coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Define function to retrieve GO term name and definition from EBI OLS
get_go_term_info <- function(go_id) {
  encoded_id <- gsub(":", "_", go_id)
  url <- paste0("https://www.ebi.ac.uk/ols/api/ontologies/go/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F", encoded_id)
  res <- GET(url)
  if (res$status_code == 200) {
    term_info <- content(res, as = "parsed", type = "application/json")
    data.frame(
      GO_ID = go_id,
      Name = term_info$label,
      Definition = if (!is.null(term_info$description)) term_info$description[[1]] else NA,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(GO_ID = go_id, Name = NA, Definition = NA)
  }
}


go_df <- do.call(rbind, lapply(go_ids, get_go_term_info))
print(go_df)
write.csv(go_df, "top10_go_annotations.csv", row.names = FALSE)
#-----------------------------------

