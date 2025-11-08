# ==============================================================================
# DEEP CLEANING & QUALITY CONTROL SCRIPT FOR COMBINED H3/N2 FILES
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP AND FILE PATHS (CRITICALLY REVISED)
# ------------------------------------------------------------------------------
setwd("/Users/joanna/conc")  # <- Ensure this is your working directory

# !!! EDIT THIS LINE for each run (H3 or N2) !!!
target_gene_file <- "Combined_H3_HA_sequences.fasta"
# target_gene_file <- "Combined_N2_NA_sequences.fasta" 

# Output file names
output_prefix <- gsub("\\.fasta$", "", target_gene_file)
output_clean <- paste0(output_prefix, "_CLEAN_FINAL.fasta")
metadata_out <- paste0(output_prefix, "_metadata_skeleton.csv")
metadata_final_map <- paste0(output_prefix, "_metadata_final_map.csv")

# QC Thresholds (REVISED TO PREVENT ZERO SEQUENCES)
# Loosen length filter to 70% to match typical HA/NA CDS lengths (~1400-1800 bp)
MIN_LENGTH_PCT <- 0.70
# Loosen N-content filter significantly to 5% to keep more sequences
MAX_PCT_N <- 0.05
# Set stop removal to FALSE initially to save all remaining data for inspection
REMOVE_INTERNAL_STOPS <- FALSE 


# ------------------------------------------------------------------------------
# 2. LOAD LIBRARIES (assuming Biostrings is installed)
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(readr) # Explicitly loaded
})


# ------------------------------------------------------------------------------
# 3. READ FASTA AND INITIAL METADATA EXTRACTION
# ------------------------------------------------------------------------------
cat("\n--- Reading Data ---\n")
if (!file.exists(target_gene_file)) {
  stop("Error: Target file not found at path: ", target_gene_file)
}
f_combined <- readDNAStringSet(target_gene_file)
cat("✅ Sequences loaded:", length(f_combined), "\n")


# Standardize headers and extract accession/date
headers <- names(f_combined)

parse_header <- function(hdr) {
  acc <- str_split(hdr, "\\s|\\|")[[1]][1]
  acc <- str_remove_all(acc, "^>")
  date_match <- str_extract(hdr, "\\d{4}-\\d{2}-\\d{2}|\\d{4}/\\d{2}/\\d{2}")
  
  tibble(
    header = hdr, 
    accession = acc, 
    date = if (is.na(date_match)) NA_character_ else gsub("/", "-", date_match)
  )
}

metadata <- headers %>% 
  lapply(parse_header) %>% 
  bind_rows() %>%
  mutate(seqname = paste0("seq", row_number())) 

names(f_combined) <- metadata$seqname
metadata$seqname_dna <- names(f_combined) 

write_csv(metadata, metadata_out)
cat("Metadata skeleton written to:", metadata_out, "\n")


# ------------------------------------------------------------------------------
# 4. REMOVE EXACT DUPLICATE SEQUENCES (DEDUPLICATION)
# ------------------------------------------------------------------------------
seq_chars <- as.character(f_combined)
dup_idx <- duplicated(seq_chars)
cat("\n--- Deduplication ---\n")
cat("Exact duplicate sequences found:", sum(dup_idx), "\n")

f_nodup <- f_combined[!dup_idx]
metadata_nodup <- metadata %>% filter(!dup_idx)

cat("After deduplication:", length(f_nodup), "sequences remaining.\n")


# ------------------------------------------------------------------------------
# 5. FILTER BY LENGTH AND AMBIGUOUS BASES (N)
# ------------------------------------------------------------------------------
cat("\n--- Length and Ambiguous Base Filtering ---\n")

lens <- width(f_nodup)
max_len <- max(lens)
min_allowed <- round(MIN_LENGTH_PCT * max_len)

cat(paste0("Max length: ", max_len, " bp. Minimum required length: ", min_allowed, " bp (", MIN_LENGTH_PCT*100, "%).\n"))

count_N <- vcountPattern("N", f_nodup, fixed=FALSE) + vcountPattern("n", f_nodup, fixed=FALSE)
pct_N <- count_N / lens

keep_len <- lens >= min_allowed
keep_N   <- pct_N <= MAX_PCT_N

# Tally results
cat("Sequences removed by length (too short):", sum(!keep_len), "\n")
cat(paste0("Sequences removed by N-content (>", MAX_PCT_N*100, "% Ns):"), sum(!keep_N), "\n")

# Apply both filters
keep_idx <- which(keep_len & keep_N)
f_filtered <- f_nodup[keep_idx]
metadata_filtered <- metadata_nodup %>% filter(seqname_dna %in% names(f_filtered))

cat("Final sequences after length/N filtering:", length(f_filtered), "\n")


# ------------------------------------------------------------------------------
# 6. TRANSLATE AND DETECT INTERNAL STOP CODONS (FRAMESHIFT/ERROR QC)
# ------------------------------------------------------------------------------
cat("\n--- Translation and Stop Codon QC ---\n")

# Translate all in frame 1 (Standard Code)
translated <- translate(f_filtered, if.fuzzy.codon="X") # Returns AAStringSet

# helper: detect internal stops (not counting a single terminal stop at the very end)
has_internal_stop_raw <- sapply(as.character(translated), function(pep){
  # find '*' positions
  stars <- which(strsplit(pep, "")[[1]] == "*")
  if (length(stars) == 0) return(FALSE)
  # allow a terminal stop only at last position
  if (length(stars) == 1 && stars == nchar(pep)) return(FALSE)
  return(TRUE)
})

# CRITICAL FIX: Ensure this is a simple logical vector for sum() and which()
has_internal_stop <- as.logical(as.vector(has_internal_stop_raw))

cat("Sequences flagged with internal stops (potential frameshifts/errors):", sum(has_internal_stop), "\n")


# --- CONDITIONAL REMOVAL BLOCK (FIXED) ---
if (REMOVE_INTERNAL_STOPS) {
  # This line now safely uses the logical vector
  keep_final_idx <- which(!has_internal_stop) 
  f_final <- f_filtered[keep_final_idx]
  # CRITICAL LINE: Defining the metadata_final object
  metadata_final <- metadata_filtered %>% filter(seqname_dna %in% names(f_final)) 
  
  cat("Final sequences kept (after stop codon removal):", length(f_final), "\n")
} else {
  # If not removing stops, the filtered set is the final set
  f_final <- f_filtered
  metadata_final <- metadata_filtered
  cat("Internal stop codon removal skipped (REMOVE_INTERNAL_STOPS=FALSE).\n")
}


# ------------------------------------------------------------------------------
# 7. RENAME HEADERS AND SAVE FINAL FILES (FIXED)
# ------------------------------------------------------------------------------
cat("\n--- Final Output ---\n")

if (length(f_final) > 0) {
  # Final Header Creation: Accession|Date
  metadata_final <- metadata_final %>% mutate(
    # If accession is NA, use the unique internal seqname
    accession_clean = ifelse(is.na(accession) | accession == "", seqname_dna, accession),
    # Replace NA dates with an empty string
    date_clean = ifelse(is.na(date), "", date)
  )
  
  new_headers <- paste0(metadata_final$accession_clean, "|", metadata_final$date_clean)
  names(f_final) <- new_headers
  
  # Save the cleaned FASTA file (ready for Alignment/Phylogeny)
  writeXStringSet(f_final, filepath = output_clean, format = "fasta")
  cat("💾 Cleaned FASTA saved to:", output_clean, " (", length(f_final), " sequences)\n")
  
  # Export the final metadata mapping
  metadata_final_out <- metadata_final %>% select(
    new_header = accession_clean, 
    date_for_tree = date_clean,
    original_full_header = header, 
    original_accession = accession,
    original_date = date,
    seqname_internal = seqname_dna # internal safe name
  )
  write_csv(metadata_final_out, metadata_final_map)
  cat("💾 Final metadata mapping written to:", metadata_final_map, "\n")
} else {
  cat("WARNING: No sequences remained after filtering. No final FASTA file written.\n")
}
cat("\n✅ QC Pipeline finished for:", target_gene_file, "\n")
