#!/usr/bin/env Rscript
# Usage: 
# 1. Read FASTA from TE library
# 2. Map to RepeatMasker .out file
# 3. Apply filters
# 4. Export passed entries as TSV for downstream analysis

# Arugment for the script in terminak: Rscript Generate_consensus.r (directory to RM.out file) (directory to TE library.fa)

# Libraries

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Biostrings)
})


# Helper functions:Clean fasta name

clean_numeric <- function(x) as.numeric(gsub("[^0-9.-]", "", x))


# Arguments

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) stop("Usage: Rscript Generate_consensus.R <rm_out_file> <fasta_file> [table_out] [fasta_out]")

rm_out_file <- args[1]
fasta_file <- args[2]
table_out <- ifelse(length(args) >= 3, args[3], "filtered_TE_table.tsv")
fasta_out <- ifelse(length(args) >= 4, args[4], "consensi.fa")
filtered_out_file <- "filtered_out_TE_table.tsv"

# Check files exist
if(!file.exists(rm_out_file)) stop(paste("RepeatMasker file not found:", rm_out_file))
if(!file.exists(fasta_file)) stop(paste("FASTA file not found:", fasta_file))


# Step 1: Read RepeatMasker .out

df <- fread(rm_out_file, skip = 3, fill = TRUE, quote = "")

# Remove completely empty columns
df <- df[, lapply(.SD, function(x) if(all(is.na(x))) NULL else x)]

# Assign known column names
known_cols <- c("SW_score","Perc_div","Perc_del","Perc_ins",
                "Scaffold","Query_start","Query_end","Position_from_end",
                "Strand","Rep_name","Rep_family","Rep_start","Rep_end",
                "Left_in_repeat","ID")

colnames(df)[1:length(known_cols)] <- known_cols

cat("Step 1: Rows in RepeatMasker file:", nrow(df), "\n")


# Step 2: Filter non-TE repeats

non_te_patterns <- c("Simple_repeat","Low_complexity","Satellite",
                     "scRNA","snRNA","srpRNA","tRNA","rRNA")

df <- df[!grepl(paste(non_te_patterns, collapse="|"), df$Rep_family), ]

cat("Step 2: Rows after removing non-TE repeats:", nrow(df), "\n")


# Step 3: Clean numeric columns

numeric_cols <- intersect(c("SW_score","Perc_div","Perc_del","Perc_ins",
                            "Query_start","Query_end","Position_from_end",
                            "Rep_start","Rep_end","Left_in_repeat"),
                          colnames(df))

df[, (numeric_cols) := lapply(.SD, clean_numeric), .SDcols = numeric_cols]


# Step 4: Fix negative strand hits

df[, c("Rep_start", "Rep_end") := .(pmin(Rep_start, Rep_end),
                                   pmax(Rep_start, Rep_end))]


# Step 5: Read FASTA

all_seqs <- readDNAStringSet(fasta_file)

fa_df <- data.frame(
  Rep_name_raw = names(all_seqs),
  sequence = as.character(all_seqs),
  stringsAsFactors = FALSE
)

# Extract Rep_name BEFORE '#'
fa_df$Rep_name <- sapply(strsplit(fa_df$Rep_name_raw, "#"), `[`, 1)

# Clean names to match RepeatMasker
fa_df$Rep_name <- trimws(fa_df$Rep_name)
fa_df$Rep_name <- gsub("[()]", "", fa_df$Rep_name)
fa_df$Rep_name <- gsub("-", "_", fa_df$Rep_name)

cat("Step 5: Number of sequences in FASTA:", nrow(fa_df), "\n")


# Step 6: Clean RepeatMasker names

df$Rep_name <- trimws(df$Rep_name)
df$Rep_name <- gsub("[()]", "", df$Rep_name)
df$Rep_name <- gsub("-", "_", df$Rep_name)


# Step 7: Merge RepeatMasker with FASTA

df <- df %>% left_join(fa_df[, c("Rep_name", "sequence")], by = "Rep_name")

cat("Step 6: Rows with matched FASTA sequences:", sum(!is.na(df$sequence)), "\n")

df <- df %>% filter(!is.na(sequence))


# Step 8: Create unified naming, pasting the TE family back

df <- df %>% mutate(
  Rep_name_final = paste(Rep_name, Rep_family, sep = "#")
)

#Clean names
df$Rep_name_final <- gsub("/", "_", df$Rep_name_final)
df$Rep_name_final <- gsub("\\s+", "_", df$Rep_name_final)


# Step 9: Compute metrics required in the filter (coverage and indel)

df <- df %>% mutate(
  consensus_length = nchar(sequence),
  alignment_length = Rep_end - Rep_start + 1,
  coverage = pmax(0, alignment_length / pmax(consensus_length, 1)),
  indel = Perc_del + Perc_ins
)


# Step 10: Apply the filter

filtered_df <- df %>% filter(
  Perc_div < 20,
  coverage > 0.8,
  indel < 10,
  SW_score > 225
)

filtered_out_df <- df %>% filter(
  !(Perc_div < 20 & coverage > 0.8 & indel < 10 & SW_score > 225)
)

cat("Rows passing filter:", nrow(filtered_df), "\n")
cat("Rows failing filter:", nrow(filtered_out_df), "\n")


# Step 11: Write tables

write.table(filtered_df, table_out, sep="\t", quote=FALSE, row.names=FALSE)
write.table(filtered_out_df, filtered_out_file, sep="\t", quote=FALSE, row.names=FALSE)


# Step 12: Deduplicate sequences (one consensus would create many copies)

unique_seqs_df <- filtered_df %>% distinct(sequence, .keep_all = TRUE)

cat("Removing redundant sequences:",
    nrow(filtered_df) - nrow(unique_seqs_df),
    "duplicates removed\n")


# Step 13: Write FASTA

seqs_unique <- DNAStringSet(unique_seqs_df$sequence)
names(seqs_unique) <- unique_seqs_df$Rep_name_final

writeXStringSet(seqs_unique, fasta_out)

cat("Done! Filtered tables written and deduplicated FASTA with",
    length(seqs_unique), "sequences.\n")