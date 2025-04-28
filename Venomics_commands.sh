
# -------------------------
# Bash Script: Transcriptome Analysis Pipeline
# -------------------------

# Step 1: Predict ORFs using TransDecoder
TransDecoder.LongOrfs -t GIOJ01.1.fsa_nt
TransDecoder.Predict -t GIOJ01.1.fsa_nt --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

# Step 2: Search Pfam domains
hmmscan --cpu 4 --domtblout pfam.domtblout Pfam-A.hmm longest_orfs.pep > pfam.log

# Step 3: BLASTp against UniProt
blastp -query longest_orfs.pep -db uniprot_sprot.fasta -outfmt 6 -evalue 1e-5 -num_threads 4 -out blastp.outfmt6

# Step 4: Run eggNOG-mapper
emapper.py -i GIOJ01.1.fsa_nt.transdecoder.pep --itype proteins -m diamond -o gioj_annotation --cpu 4
