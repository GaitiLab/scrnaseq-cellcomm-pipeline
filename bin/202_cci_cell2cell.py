#!/usr/bin/env python3

# Import modules
import argparse
import logging
import multiprocessing
import os
import pickle
import time
from argparse import ArgumentParser as AP
from os.path import abspath
from pathlib import Path

import cell2cell as c2c
import pandas as pd
import scprep


def get_args():
    # Script description
    description = """Infer CCIs with Cell2Cell"""

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
        help="Path to custom database with interactions (csv)",
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

    parser.add_argument("--version", action="version", version="0.1.0")
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


def parse_path(p):
    if os.path.islink(p):
        return os.readlink(p)
    else:
        return p


def run_cell2cell(
    sample_id: str,
    input_dir: str,
    meta_path: str,
    interactions_db: str,
    annot: str,
    output_dir: str,
    n_perm: int = 1000,
):
    logging.info(f"Current sample id: {sample_id}")
    logging.info("Load RNAseq data (10x format)...")
    if os.path.exists(args.input_dir):
        rnaseq_df = scprep.io.load_10X(input_dir).T
    else:
        raise Exception("Not a valid path for the input file")
    # Set genes as index
    rnaseq_df = rnaseq_df.rename_axis("index", axis=1)

    logging.info("Load database with CCIs...")
    lr_pairs = pd.read_csv(interactions_db)
    lr_pairs = lr_pairs.astype(str)

    # Metadata for the single cells
    logging.info("Load metadata...")
    meta_df = pd.read_csv(meta_path, index_col=0)
    meta_df.index.name = "index"
    meta_df = meta_df.loc[meta_df.index.isin(rnaseq_df.columns)]

    # Cell-cell Interactions and Communication Analysis
    # The pipeline integrates the RNA-seq and PPI datasets by using the analysis setups.
    # It generates an interaction space containing an instance for each sample/cell type, containing the values assigned to each protein in the PPI list given the setups for computing the CCI and CCC scores.

    logging.info("Infer CCIs")
    interactions = c2c.analysis.SingleCellInteractions(
        rnaseq_data=rnaseq_df,
        ppi_data=lr_pairs,
        metadata=meta_df,
        interaction_columns=("source_genesymbol", "target_genesymbol"),
        communication_score="expression_gmean",
        cci_score="bray_curtis",
        cci_type="directed",
        aggregation_method="average",
        barcode_col="index",
        celltype_col=annot,
        complex_sep="_",
        verbose=True,
    )

    # **Compute communication scores for each PPI or LR pair**
    interactions.compute_pairwise_communication_scores()
    logging.info("Perform permutation analysis...")
    interactions.permute_cell_labels(
        evaluation="communication",
        permutations=n_perm,
        fdr_correction=True,
        verbose=True,
    )
    logging.info("Save results...")
    # Save p-values
    interactions.ccc_permutation_pvalues.to_csv(
        f"{output_dir}/cell2cell__{sample_id}__pvalues.csv", sep="\t"
    )

    # Save interaction scores
    interactions.interaction_space.interaction_elements["communication_matrix"].to_csv(
        f"{output_dir}/cell2cell__{sample_id}__interaction_scores.csv",
        sep="\t",
    )

    logging.info("Save interactions object as pickle...")
    with open(f"{output_dir}/cell2cell__{sample_id}.pickle", "wb") as handle:
        pickle.dump(interactions, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logging.info("COMPLETED")


def main(args):
    # Setup logging
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(message)s")

    run_cell2cell(
        sample_id=args.sample_id,
        input_dir=args.input_dir,
        meta_path=args.meta,
        annot=args.annot,
        output_dir=args.output_dir,
        n_perm=args.n_perm,
        interactions_db=args.interactions_db,
    )


if __name__ == "__main__":
    args = get_args()
    st = time.time()
    main(args)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")
