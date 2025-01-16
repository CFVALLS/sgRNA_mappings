#!/bin/bash

#function to display progress messages
progress_msg() {
    echo "$(date +"%Y-%m-%d %H:%M:%S") - $1"
}

# Function to display usage
usage() {
    echo "Usage: $0 -g <genome_index> -s <sgRNA_library.fa> -a <annotation.gtf> -o <output_dir> [-p <output_prefix>]"
    echo "  -g: Path to the Bowtie2 genome index" # /data/shared/assemblies/GRCh38_110/bowtie2_index/GRCH38_110_index
    echo "  -s: Path to the sgRNA library FASTA file"
    echo "  -a: Path to the annotation GTF file"
    echo "  -o: Output directory"
    echo "  -p: Output prefix (optional, default: 'output')"
    exit 1
}

# Default values
OUTPUT_PREFIX="output"
CLEAN_AFTER=false

# Parse command line arguments
while getopts ":g:s:a:o:p:" opt; do
    case $opt in
        g) GENOME_INDEX="$OPTARG";;
        s) SGRNA_LIBRARY="$OPTARG";;
        a) ANNOTATION_GTF="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        p) OUTPUT_PREFIX="$OPTARG";;
        \?) echo "Invalid option -$OPTARG" >&2; usage;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage;;
    esac
done

# Check if required arguments are provided
if [ -z "$GENOME_INDEX" ] || [ -z "$SGRNA_LIBRARY" ] || [ -z "$ANNOTATION_GTF" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Mapping sgRNAs to the genome
progress_msg "Starting Bowtie2 mapping of sgRNAs to the genome..."

bowtie2 -x "$GENOME_INDEX" -f "$SGRNA_LIBRARY" -k 100 \
    -N 0 --gbar 1 --ma 2 --mp 0,0 --rdg 999,999 --rfg 999,999 --score-min "C,2,0" \
    -S "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping.sam" \
    --un "${OUTPUT_DIR}/${OUTPUT_PREFIX}_unmapped.fa" \
    --threads $(nproc)
progress_msg "Bowtie2 mapping completed."

# Convert SAM to BAM, sort, and index
progress_msg "Converting SAM to BAM, sorting, and indexing..."
samtools view -bS "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping.sam" | \
    samtools sort -@ $(nproc) -o "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping_sorted.bam"
samtools index "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping_sorted.bam"
progress_msg "BAM conversion, sorting, and indexing completed."

# Convert BAM to BED
progress_msg "Converting BAM to BED..."
bedtools bamtobed -i "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping_sorted.bam" > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping.bed"
progress_msg "BAM to BED conversion completed."

# Annotate mappings using BEDTools
# -wa 	Write the original entry in A for each overlap.
# -wb 	Write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r.
# loj left-out-join 
progress_msg "Annotating mappings with BEDTools..."
bedtools intersect -a "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping.bed" -b "$ANNOTATION_GTF" -wa -wb -loj > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_annotated_mappings.bed"
progress_msg "Annotation completed."

# Count multi-mapping reads
progress_msg "Counting multi-mapping reads..."
cut -f4 "${OUTPUT_DIR}/${OUTPUT_PREFIX}_annotated_mappings.bed" | sort | uniq -c > "${OUTPUT_DIR}/${OUTPUT_PREFIX}_mapping_counts.txt"
progress_msg "Multi-mapping read count completed."

# Clean up intermediate files
if [[ "$CLEAN_AFTER" == true ]]; then  
    progress_msg "Cleaning up intermediate files..."
    rm "${OUTPUT_DIR}/${OUTPUT_PREFIX}_genome_mapping.sam"
    progress_msg "Cleanup completed."
fi

progress_msg "All processes completed successfully."
echo "Output files in ${OUTPUT_DIR}:"
echo "  - ${OUTPUT_PREFIX}_genome_mapping_sorted.bam"
echo "  - ${OUTPUT_PREFIX}_genome_mapping_sorted.bam.bai"
echo "  - ${OUTPUT_PREFIX}_genome_mapping.bed"
echo "  - ${OUTPUT_PREFIX}_annotated_mappings.bed"
echo "  - ${OUTPUT_PREFIX}_unmapped.fa"
echo "  - ${OUTPUT_PREFIX}_mapping_counts.txt"
