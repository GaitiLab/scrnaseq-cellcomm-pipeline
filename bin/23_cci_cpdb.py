#!/usr/bin/env python3
import argparse
import logging
import multiprocessing
import os
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import pandas as pd
import scanpy as sc
import scprep
import utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method


def get_args():
    # Script description
    description = """Infer CCIs with CellPhoneDB v5"""

    # Add parser
    parser = AP(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-f", "--input_dir", type=str, help="Path to 10X data directory"
    )
    parser.add_argument(
        "-p",
        "--n_perm",
        type=int,
        default=1000,
        help="Number of permutations for permutation testing",
    )
    parser.add_argument(
        "-db",
        "--interactions_db",
        type=str,
        help="Path to custom database with interactions (zip)",
    )
    parser.add_argument(
        "-a",
        "--annot",
        type=str,
        help="Column in metadata containing the cell type labels",
    )
    parser.add_argument("-id", "--sample_id", type=str, default="")
    parser.add_argument(
        "-o", "--output_dir", type=str, default="output", help="Output directory"
    )
    parser.add_argument("-m", "--meta", type=str, help="Path to metadata file (CSV)")
    parser.add_argument(
        "-nc",
        "--n_cores",
        type=int,
        default=4,
        help="Number of cores to use for parallelization",
    )
    parser.add_argument(
        "-mp",
        "--min_pct",
        type=float,
        default=0.1,
        help="Minimum percentage of cells expressing a gene [0.0-1.0], a.k.a. proportion",
    )

    parser.add_argument(
        "--nf-process-id",
        type=str,
        help="Nextflow process ID",
        default=None,
        dest="nf_process_id",
    )
    arg = parser.parse_args()
    arg.output_dir = abspath(arg.output_dir)
    arg.input_dir = abspath(arg.input_dir)
    arg.interactions_db = abspath(arg.interactions_db)

    if (arg.output_dir != "") & (not os.path.isdir(arg.output_dir)):
        # Create an empty folder for TF records if folder doesn't exist
        # arg.output_dir = Path(arg.output_dir)
        Path(arg.output_dir).parent.mkdir(parents=True, exist_ok=True)

        if arg.n_cores is None:
            arg.n_cores = multiprocessing.cpu_count()
    return arg


def run_cpdb(
    sample_id: str,
    input_dir: str,
    meta_path: str,
    interactions_db: str,
    annot: str,
    output_dir: str,
    n_perm: int = 1000,
    min_pct: float = 0.1,
    n_cores: int = 4,
):
    logging.info(f"Current sample id: {sample_id}")

    logging.info("Load metadata...")
    metadata = pd.read_csv(meta_path, sep=",", index_col=0)

    logging.info("Load RNAseq data (10x format)...")
    rnaseq = scprep.io.load_10X(input_dir)
    rnaseq.index.name = "barcode_sample"
    rnaseq.columns.name = "gene_ids"

    logging.info("Create AnnData object for CellPhoneDB...")
    adata = sc.AnnData(rnaseq, rnaseq.index.to_frame(), rnaseq.columns.to_frame())

    logging.info("Subsetting metadata to match the rnaseq data...")
    metadata = metadata.loc[adata.obs.index]

    logging.info(
        "Renaming the cell type column to match the CellPhoneDB format 'cell_type'..."
    )
    metadata = metadata.rename(columns={annot: "cell_type"})

    logging.info("Writing AnnData object to disk...")
    adata.write_h5ad(f"{output_dir}/{sample_id}_counts.h5ad")

    logging.info("Writing metadata to disk...")
    metadata.to_csv(f"{output_dir}/{sample_id}_metadata.tsv", sep="\t", index=True)

    logging.info("Setup paths...")
    counts_file_path = f"{output_dir}/{sample_id}_counts.h5ad"
    meta_file_path = f"{output_dir}/{sample_id}_metadata.tsv"

    logging.info("Run CellPhoneDB...")
    cpdb_results = cpdb_statistical_analysis_method.call(
        # mandatory: CellphoneDB database zip file.
        cpdb_file_path=interactions_db,
        # mandatory: tsv file defining barcodes to cell label.
        meta_file_path=meta_file_path,
        # mandatory: normalized count matrix.
        counts_file_path=counts_file_path,
        # defines the gene annotation in counts matrix.
        counts_data="hgnc_symbol",
        # active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
        # microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
        # optional: whether to score interactions or not.
        score_interactions=True,
        # denotes the number of shufflings performed in the analysis.
        iterations=n_perm,
        # defines the min % of cells expressing a gene for this to be employed in the analysis.
        threshold=min_pct,
        # number of threads to use in the analysis.
        threads=n_cores,
        # debug randome seed. To disable >=0.
        # debug_seed=42,
        # Sets the rounding for the mean values in significan_means.
        result_precision=7,
        # P-value threshold to employ for significance.
        # Keep all interactions, filtering a posteriori
        pvalue=1.1,
        # To enable subsampling the data (geometri sketching).
        subsampling=False,
        # subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
        # Number of componets to subsample via geometric skectching (dafault: 100).
        # subsampling_num_pc=100,
        # subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
        # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        separator="|",
        # Saves all intermediate tables employed during the analysis in pkl format.
        debug=False,
        # Path to save results.
        output_path=output_dir,
        # Replaces the timestamp in the output files by a user defined string in the  (default: None).
        output_suffix=f"_{sample_id}",
    )
    logging.info("COMPLETED")


def main(args):
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    run_cpdb(
        sample_id=args.sample_id,
        input_dir=args.input_dir,
        meta_path=args.meta,
        interactions_db=args.interactions_db,
        annot=args.annot,
        output_dir=args.output_dir,
        n_perm=args.n_perm,
        min_pct=args.min_pct,
        n_cores=args.n_cores,
    )

    if args.nf_process_id is not None:
        utils.generate_versions_yml(
            ["scprep", "scprep", "scanpy", "pandas", "cellphonedb"],
            task_id=args.nf_process_id,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
