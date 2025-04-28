# Set column names based on BLAST outfmt 6
colnames <- c("query_id", "subject_id", "perc_identity", "alignment_length",
              "mismatches", "gap_opens", "q_start", "q_end", 
              "s_start", "s_end", "evalue", "bit_score")

# Load your BLAST results
blastp <- read.table("C:/Users/menna/OneDrive/Documents/IU LUDDY BIOINFORMATICS/Comp analysis high TP biomedical data/Project/Results/TransDecoder/venom_vs_transdecoder.blastp", header = FALSE, sep = "\t", col.names = colnames)

# Keep the best hit per venom protein (lowest evalue, highest bit score)
library(dplyr)
best_hits <- blastp %>%
  group_by(query_id) %>%
  arrange(evalue, desc(bit_score)) %>%
  slice(1) %>%
  ungroup()

# Save to CSV
write.csv(best_hits, "best_blastp_hits.csv", row.names = FALSE)
