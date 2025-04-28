# HTP_Project_Venomics
Functional annotation and transcriptomic validation of venom proteins from Acanthoscurria rondoniae using a multiomics computational pipeline. Includes coding sequence prediction, GO/KEGG/COG annotation, venom protein matching, and enrichment analysis.

## Script Descriptions
This repository contains the core scripts used for the computational analysis and functional annotation of the Acanthoscurria rondoniae venom transcriptome.

---

### 1. Venomics_commands.sh
Bash script used for:
- ORF prediction with TransDecoder
- Domain identification with Pfam (hmmscan)
- Sequence similarity search with BLASTp
- Functional annotation using eggNOG-mapper

**Required input:**
- `GIOJ01.1.fsa_nt` (transcriptome FASTA)

---

### 2. best_blastp_hits.R
Filters BLASTp results to retain best transcript matches for each curated venom protein.

**Input:**
- `blastp_results.csv`  
**Output:**
- `best_blastp_hits.csv`

---

### 3. GO and KEGG enrichment_eggNOG_result.R
R script for functional profiling of predicted proteins:
- Parses eggNOG annotations
- Extracts and visualizes top GO terms and KEGG orthologs
- Maps KEGG KOs to pathway names using KEGGREST

**Input:**
- `gioj_annotation.emapper.annotations`

---

### 4. project_script.R
Parses mzIdentML output and matches peptides to curated venom proteins.

**Input:**
- `peptides_1_1_0_Ar_proteomics.mzid.xml`  
- `venom_table2_curated_proteins.csv`  
**Output:**
- `matched_peptides_to_table2.csv`

---

### 5. GO_Term_Protein_netwrok.R
Generates a bipartite GO–Protein network:
- Visualizes relationships between GO terms and curated venom proteins
- Uses `ggraph` and `igraph`

**Input:**
- `gene_list_for_topGO.csv`
- `GO_enrichment_results.csv`
- `protein_to_go_mapping.csv`

**Output:**
- GO–Protein network plot (PNG)

---

## Software Dependencies

- TransDecoder v5.5.0
- BLAST+ v2.14.1
- HMMER (hmmscan) with Pfam-A.hmm
- eggNOG-mapper v2.1.9
- R (topGO, tidyverse, ggraph, KEGGREST)
- Python (pandas, seaborn, matplotlib) [original only]

---

## Reproducibility Notes

Ensure required input files are present in the working directory and modify paths in each script accordingly. All intermediate steps are tracked and outputs saved for reproducibility.
