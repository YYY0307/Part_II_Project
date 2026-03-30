#!/bin/bash


# FULL TE PIPELINE: Automated scrip to mask and Generate the alignments file (in stks.) based on RepeatMasker.out file, TE library and genome
# Dependencies: Make sure Generated_consensus.R is installed in the computer and kept in the same directory
# Usage:
# This version is created for users that installed RepeatMasker and ReportModeler through dockker/dfam-te:tools. For Rosetta users please modify directory to these applications if appropriate.
# ./full_TE_pipeline.sh <RM_out> <library_fasta> <genome_fasta> <species_name> <output_prefix>

set -e  # stop on error

RM_OUT=$1
LIB_FA=$2
GENOME_FA=$3
SPECIES_NAME=$4
PREFIX=$5

THREADS=8


# Directories

BASE_DIR="${PREFIX}_TE_pipeline"
CONSENSUS_DIR="$BASE_DIR/consensus"
RM_DIR="$BASE_DIR/RepeatMasker"
STKS_DIR="$BASE_DIR/stks"

mkdir -p "$CONSENSUS_DIR"
mkdir -p "$RM_DIR"
mkdir -p "$STKS_DIR"


echo "Running TE pipeline for: $SPECIES_NAME"


# Step 1: Run R consensus pipeline

echo "Step 1: Generating filtered consensus library"

Rscript Generate_consensus.R \
  "$RM_OUT" \
  "$LIB_FA" \
  "$CONSENSUS_DIR/filtered_TE_table.tsv" \
  "$CONSENSUS_DIR/consensi.fa"

CONSENSI_FA="$CONSENSUS_DIR/consensi.fa"


# Step 2: Convert genome to 2bit

echo "Step 2: Converting genome to 2bit"

docker run -it --rm \
  -v "$PWD":/data \
  -w /data \
  dfam/tetools:latest \
  faToTwoBit "$GENOME_FA" "${GENOME_FA%.fa}.2bit"

TWO_BIT="${GENOME_FA%.fa}.2bit"


# Step 3: Run RepeatMasker to mask genomes

echo "Step 3: Running RepeatMasker"

RepeatMasker \
  -pa "$THREADS" \
  -lib "$CONSENSI_FA" \
  -a \
  -dir "$RM_DIR" \
  "$GENOME_FA"

ALIGN_FILE="$RM_DIR/$(basename ${GENOME_FA}).align"


# Step 4: Generate seed alignments

echo "Step 4: Running generateSeedAlignments"

docker run -it --rm \
  -v "$PWD":/data \
  -w /data/"$STKS_DIR" \
  dfam/tetools:latest \
  generateSeedAlignments.pl \
    -taxon "$SPECIES_NAME" \
    -assemblyFile /data/"$TWO_BIT" \
    /data/"$ALIGN_FILE" \
    > generateSeedAlignments.log 2>&1


# Step 5: Run Linup on .stk files

echo "Step 5: Running Linup"

docker run -it --rm \
  -v "$PWD":/data \
  -w /data/"$STKS_DIR" \
  dfam/tetools:latest \
  bash -c 'for infile in *.stk; do
    Linup "$infile" -consensus > "${infile}_con.fa";
    Linup "$infile" -fasta     > "${infile}_elements.fa";
    echo "done: $infile";
  done'

echo "✅ PIPELINE COMPLETE for $SPECIES_NAME"
echo "Outputs in: $BASE_DIR"
