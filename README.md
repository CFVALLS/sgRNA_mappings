
## Overview - So i dont forget
This pipeline consists of a series of scripts and tools designed to map sgRNA libraries to a reference genome, identify overlaps with annotated ORFs (Open Reading Frames), and aggregate the results for downstream analysis. It is particularly useful for analyzing CRISPR knockout screens.

### Scripts in the Pipeline
1. **`map_sgRNA_Lib.sh`**: Maps sgRNA sequences to a reference genome using Bowtie2 and annotates the mappings.
2. **`map_CRISPRko_guides.sh`**: Facilitates the processing of sgRNA mappings, overlaps with ORFs, and filtering.
3. **`aggregate_mappings.py`**: Aggregates the sgRNA mappings, gene annotations, and ORF overlaps into a summary CSV file.

---

## Inputs

### Required Files
1. **Genome Index** (`GENOME_INDEX`): Bowtie2 index of the reference genome (e.g., `GRCH38_110_index`).
2. **sgRNA FASTA File** (`SGRNA_FASTA`): FASTA file containing sgRNA sequences and associated gene annotations (e.g., `CRISPRko_Brunello_guide_gene_map.fasta`).
3. **Annotation GTF File** (`ANNOTATION_GTF`): GTF file with gene and transcript annotations (e.g., `gencode.v47.Riboseq_ORFs.complete.sorted.gtf`).
4. **Ribo-seq ORFs BED File** (`RIBOSEQ_ORFS`): BED file containing annotated ribosomal ORFs (e.g., `Ribo-seq_ORFs_nochar_tab.bed`).

### Optional Parameters
- **Output Directory** (`OUTPUT_DIR`): Directory to store intermediate and final results (default is predefined).
- **Prefix** (`PREFIX`): Prefix for naming intermediate files (default: `crispr_screen`).

---

## Outputs

### Generated Files
1. **Genome Mapping BED File** (`crispr_screen_genome_mapping.bed`): Mapped sgRNA positions in the genome.
2. **ORF Overlaps File** (`sgRNA_orf_overlaps_fixed.bed`): Overlap data between sgRNA mappings and Ribo-seq ORFs.
3. **Filtered Overlaps File** (`sgRNA_orf_overlaps_filtered.bed`): Overlaps with positive values filtered for relevance.
4. **Aggregated Summary CSV** (`mapping_summary.csv`): Summary table combining sgRNA mappings, gene annotations, and ORF overlaps.

---

## Workflow Steps

1. **Map sgRNA Library**:
   The `map_sgRNA_Lib.sh` script aligns sgRNA sequences from the FASTA file to the reference genome using Bowtie2. The output is a BED file containing genome mappings.

2. **Identify ORF Overlaps**:
   Using `bedtools intersect`, the mapped sgRNA positions are compared against the Ribo-seq ORF annotations to find overlaps. Results are stored in `sgRNA_orf_overlaps_fixed.bed`.

3. **Filter Overlaps**:
   Using `awk`, overlaps with positive values (indicating actual intersections) are filtered and saved to `sgRNA_orf_overlaps_filtered.bed`.

4. **Aggregate Mappings**:
   The `aggregate_mappings.py` script combines the sgRNA FASTA file, genome mapping file, gene annotations, and filtered ORF overlaps into a single summary CSV file for further analysis.

---

## Example Usage

### Bash Script
Run the pipeline with the following command:
```bash
bash map_CRISPRko_guides.sh
