#!/bin/bash

# Set variables for file paths
GENOME_INDEX=~/shared/assemblies/GRCh38_110/bowtie2_index/GRCH38_110_index
SGRNA_FASTA=/data/cfvall/projects/reprocessing_Sanson_et_al_20250801/data/processed_files/CRISPRko_Brunello_guide_gene_map.fasta
ANNOTATION_GTF=/data/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/gencode.v47.Riboseq_ORFs.complete.sorted.gtf
OUTPUT_DIR=/data/cfvall/projects/reprocessing_Sanson_et_al_20250801/data/mapping_sgRNA_Brunello
BEDTOOLS_INPUT=${OUTPUT_DIR}/crispr_screen_genome_mapping.bed
RIBOSEQ_ORFS=~/shared/assemblies/GENCODE_Phase1_Riboseq_ORFs/Ribo-seq_ORFs_nochar_tab.bed
BEDTOOLS_OUTPUT=${OUTPUT_DIR}/sgRNA_orf_overlaps_fixed.bed
FILTERED_OUTPUT=${OUTPUT_DIR}/sgRNA_orf_overlaps_filtered.bed
PREFIX=crispr_screen

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Step 1: Map sgRNA library using custom script
./map_sgRNA_Lib.sh \
    -g "$GENOME_INDEX" \
    -s "$SGRNA_FASTA" \
    -a "$ANNOTATION_GTF" \
    -o "$OUTPUT_DIR" \
    -p "$PREFIX"

# Step 2: Find overlaps between sgRNA mappings and Ribo-seq ORFs
bedtools intersect \
    -a "$BEDTOOLS_INPUT" \
    -b "$RIBOSEQ_ORFS" \
    -wao > "$BEDTOOLS_OUTPUT"

# Step 3: Filter results to include only overlaps with positive values
awk '$NF > 0' "$BEDTOOLS_OUTPUT" > "$FILTERED_OUTPUT"

# Print a completion message
echo "sgRNA mapping and filtering completed. Results saved to:"
echo "$FILTERED_OUTPUT"

# Step 4: Aggregate mappings
python3 aggregate_mappings.py \
    --sgrna-file "$SGRNA_FASTA" \
    --genome-mapping "$BEDTOOLS_INPUT" \
    --gene-mapping "${OUTPUT_DIR}/${PREFIX}_annotated_mappings.bed" \
    --orf-mapping "$FILTERED_OUTPUT" \
    --output "$OUTPUT_DIR/mapping_summary.csv"
