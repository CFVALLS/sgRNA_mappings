import pandas as pd
import argparse
from pathlib import Path
import logging


def setup_logging():
    """Configure logging for the script."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def read_sgrna_data(fasta_path: str) -> pd.DataFrame:
    """Read and process sgRNA data from a FASTA file."""
    logging.info("Reading sgRNA data...")
    sequences = []
    with open(fasta_path, "r") as f:
        header = None
        sequence = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append({"guide_id": header, "sgrna": sequence})
                header = line[1:]  # Remove '>' character
                sequence = None
            else:
                sequence = line
        if header is not None:  # Add the last sequence
            sequences.append({"guide_id": header, "sgrna": sequence})

    df_fasta = pd.DataFrame(sequences)
    if "::" not in df_fasta["guide_id"].iloc[0]:
        raise ValueError("The 'guide_id' field must contain 'gene::id' format.")
    df_fasta["id"] = df_fasta["guide_id"].str.split("::").str[1]
    df_fasta["gene"] = df_fasta["guide_id"].str.split("::").str[0]

    return df_fasta


def validate_file(file_path: str):
    """Ensure the file exists and is readable."""
    if not Path(file_path).exists():
        raise FileNotFoundError(f"File not found: {file_path}")


def process_genome_mappings(file_path: str) -> pd.DataFrame:
    """Process genome mapping data."""
    logging.info("Processing genome mappings...")
    validate_file(file_path)
    mappings = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "guide_id", "score", "strand"],
    )
    return mappings[mappings["score"] > 0]


def process_gene_mappings(file_path: str) -> pd.DataFrame:
    """Process gene mapping data."""
    logging.info("Processing gene mappings...")
    column_names = [
        "chr", "start", "end", "guide_id", "MAPQ", "strand",
        "chr_annotation", "source_annotation", "feature_annotation",
        "start_annotation", "end_annotation", "frame_annotation",
        "strand_annotation", "score_annotation", "additional_info",
    ]
    gene_mapping = pd.read_csv(file_path, sep="\t", header=None, names=column_names)
    gene_mapping = gene_mapping[
        (gene_mapping["MAPQ"] > 0) & (gene_mapping["feature_annotation"] == "gene")
    ]

    gene_mapping = gene_mapping.assign(
        ensembl_gene_id_annotation=gene_mapping["additional_info"].str.extract(r'gene_id "([^"]+)"'),
        ensembl_transcript_id_annotation=gene_mapping["additional_info"].str.extract(r'transcript_id "([^"]+)"'),
        gene_type=gene_mapping["additional_info"].str.extract(r'gene_type "([^"]+)"'),
        gene_name_annotation=gene_mapping["additional_info"].str.extract(r'gene_name "([^"]+)"'),
        ncORF_name_annotation=gene_mapping["additional_info"].str.extract(r'orf_id "([^"]+)"'),
    )

    unique_genes = gene_mapping.groupby("guide_id")["ensembl_gene_id_annotation"].unique().reset_index()
    unique_genes = unique_genes.rename(columns={"ensembl_gene_id_annotation": "unique_genes_targeted"})

    count_genes = gene_mapping.groupby("guide_id")["ensembl_gene_id_annotation"].nunique().reset_index()
    count_genes = count_genes.rename(columns={"ensembl_gene_id_annotation": "count_gene_targeted"})

    return gene_mapping, unique_genes, count_genes


def process_orf_data(file_path: str) -> pd.DataFrame:
    """Process ORF intersection data."""
    logging.info("Processing ORF data...")
    orfeome_intersect = pd.read_csv(file_path, sep="\t", header=None)
    orfeome_intersect = orfeome_intersect[orfeome_intersect[4] > 0]

    orf_counts = orfeome_intersect.groupby(3)[9].nunique().reset_index()
    orf_counts = orf_counts.rename(columns={3: "guide_id", 9: "count_orf_targeted"})

    orf_lists = orfeome_intersect.groupby(3)[9].agg(list).reset_index()
    orf_lists = orf_lists.rename(columns={3: "guide_id", 9: "orf_ids_targeted"})

    return orf_counts, orf_lists


def process_crispr_mappings(
    sgrna_file: str,
    genome_mapping_file: str,
    gene_mapping_file: str,
    orf_mapping_file: str,
) -> pd.DataFrame:
    """Process CRISPR screen mapping data."""
    logging.info("Processing CRISPR mappings...")
    df_sgRNAs = read_sgrna_data(sgrna_file)
    genome_mappings = process_genome_mappings(genome_mapping_file)

    df_merged = pd.merge(genome_mappings, df_sgRNAs, on="guide_id", how="left")

    mapping_counts = df_merged.groupby("guide_id")["score"].count().reset_index()
    mapping_counts = mapping_counts.rename(columns={"score": "high_quality_genomic_mappings"})
    df_merged = pd.merge(df_merged, mapping_counts, on="guide_id", how="left")

    _, unique_genes, gene_counts = process_gene_mappings(gene_mapping_file)
    df_merged = pd.merge(df_merged, gene_counts, on="guide_id", how="left")
    df_merged = pd.merge(df_merged, unique_genes, on="guide_id", how="left")

    orf_counts, orf_lists = process_orf_data(orf_mapping_file)
    df_merged = pd.merge(df_merged, orf_counts, on="guide_id", how="left")
    df_merged = pd.merge(df_merged, orf_lists, on="guide_id", how="left")

    df_merged["count_orf_targeted"] = df_merged["count_orf_targeted"].fillna(0)
    df_merged["orf_ids_targeted"] = df_merged["orf_ids_targeted"].apply(
        lambda x: x if isinstance(x, list) else []
    )

    logging.info("Processing completed successfully.")
    return df_merged


def main():
    """Main function to execute the script."""
    parser = argparse.ArgumentParser(description="Process CRISPR screen mapping data")
    parser.add_argument("--sgrna-file", required=True, help="Path to sgRNA FASTA file")
    parser.add_argument("--genome-mapping", required=True, help="Path to genome mapping file")
    parser.add_argument("--gene-mapping", required=True, help="Path to gene mapping file")
    parser.add_argument("--orf-mapping", required=True, help="Path to ORF mapping file")
    parser.add_argument("--output", required=True, help="Path to output file")

    args = parser.parse_args()
    setup_logging()

    try:
        logging.info("Starting CRISPR mapping process...")
        result_df = process_crispr_mappings(
            args.sgrna_file,
            args.genome_mapping,
            args.gene_mapping,
            args.orf_mapping,
        )
        result_df.drop(columns = ["chr","start","end", "score", "strand"], inplace = True)
        result_df.to_csv(args.output, index=False)
        logging.info(f"Results saved to {args.output}")
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise


if __name__ == "__main__":
    main()
