# ----------------------------------------
# 📦 Install Required Packages
# ----------------------------------------
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")

# ----------------------------------------
# 📚 Load Libraries
# ----------------------------------------
library(tidyverse)
library(topGO)
library(dplyr)

# ----------------------------------------
# 📂 Load Cleaned InterProScan Annotations
# ----------------------------------------
df <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/interproscan_annotations_clean.csv")

# ----------------------------------------
# 🧹 Keep Only ProteinID and GO_Terms
# ----------------------------------------
df_go <- df %>%
  dplyr::select(ProteinID, GO_Terms) %>%
  dplyr::filter(!is.na(GO_Terms)) %>%
  separate_rows(GO_Terms, sep = ";") %>%
  distinct()

# ----------------------------------------
# 🔁 Create gene-to-GO mapping
# ----------------------------------------
geneID2GO <- split(df_go$GO_Terms, df_go$ProteinID)

# ----------------------------------------
# 🧬 All Genes with GO annotations
# ----------------------------------------
allGenes <- names(geneID2GO)

# ----------------------------------------
# 🎯 Load Curated Gene List (from .csv file or manually)
# Make sure curated_genes is a CHARACTER vector
# ----------------------------------------
curated_df <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/gene_list_for_topGO.csv")
curated_genes <- curated_df$gene  # Ensure this is the correct column name

# ----------------------------------------
# 🔎 Define Interesting Genes
# ----------------------------------------
interestingGenes <- curated_genes

# ----------------------------------------
# 🎯 Create geneList Vector
# (factor with 2 levels: 0 = not interesting, 1 = interesting)
# ----------------------------------------
geneList <- factor(as.integer(allGenes %in% interestingGenes))
names(geneList) <- allGenes
levels(geneList) <- c("0", "1")

# ----------------------------------------
# 🧠 Create topGOdata Object
# ----------------------------------------
GOdata <- new("topGOdata",
              ontology = "BP",  # Use "MF" or "CC" for other ontologies
              allGenes = geneList,
              geneSel = function(x) x == 1,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO,
              description = "GO Enrichment of curated venom proteins")

# ----------------------------------------
# 🧪 Run Enrichment Test
# ----------------------------------------
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# ----------------------------------------
# 📊 Display Top GO Results
# ----------------------------------------
go_results <- GenTable(GOdata, 
                       classicFisher = resultFisher, 
                       topNodes = 15)

# View results
print(go_results)

# Optional: Save to CSV
write_csv(go_results, "GO_enrichment_results.csv")

# ----------------------------------------
# Load the topGO results
go_df <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/TopGO/GO_enrichment_results.csv")

# Select top 10 enriched terms by p-value
top_terms <- go_df %>%
  arrange(classicFisher) %>%
  slice(1:10)

# Plot
ggplot(top_terms, aes(x = reorder(Term, -classicFisher), y = -log10(classicFisher))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched GO Terms (topGO)",
       x = "GO Term", y = "-log10(p-value)") +
  theme_minimal()
