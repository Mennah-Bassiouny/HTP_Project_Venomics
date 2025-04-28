# Load required libraries
library(tidyverse)
library(ggraph)
library(igraph)
library(readr)

# Step 0: Load data
enriched <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/TopGO/GO_enrichment_results.csv")
mapping <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/protein_to_go_mapping.csv")
curated <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/gene_list_for_topGO.csv")
curated_proteins_df <- read_csv("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Datasets/venom_table2_curated_proteins.csv")

# Step 1: Prepare vectors
enriched_go_ids <- as.character(enriched$GO.ID)
curated_gene_ids <- as.character(curated$gene)

# Step 2: Clean mapping file
colnames(mapping) <- c("ProteinID", "GO")

# Step 3: Filter for enriched GO terms and curated proteins
filtered_mapping <- mapping %>%
  filter(GO %in% enriched_go_ids, ProteinID %in% curated_gene_ids)

# Step 4: Build edge list
edges <- filtered_mapping %>%
  rename(from = GO, to = ProteinID)

# Step 5: Create node table
nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(type = ifelse(name %in% enriched_go_ids, "GO Term", "Protein"))

# Step 6: Prepare label mappings
go_term_map <- enriched %>% 
  select(name = GO.ID, Label = Term)

protein_name_map <- curated_proteins_df %>%
  select(name = Accession, Label = Description)

# Step 7: Merge labels
nodes <- nodes %>%
  left_join(go_term_map, by = "name") %>%
  left_join(protein_name_map, by = "name") %>%
  mutate(Label = coalesce(Label.x, Label.y, name)) %>%
  select(name, type, Label)

# Step 8: Create graph
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Step 9: Plot and save image with enhanced visibility
p <- ggraph(graph, layout = "fr") + 
  geom_edge_link(alpha = 0.4, color = "gray50", width = 0.3) +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = Label), repel = TRUE, size = 4, color = "black") +
  scale_color_manual(values = c("GO Term" = "#1f78b4", "Protein" = "#33a02c")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "GO Term – Protein Network (topGO)",
    color = "Node Type"
  )

# Save high-resolution plot
ggsave("GO_term_protein_network_labeled_fixed.png", plot = p, width = 14, height = 12, dpi = 300)
