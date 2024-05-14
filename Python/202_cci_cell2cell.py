# Import modules
import argparse
import logging
import cell2cell as c2c
import pandas as pd
import scprep
import os
import numpy as np
import pyreadr as py
from pathlib import Path
import pickle


def parse_path(p):
    if os.path.islink(p):
        return os.readlink(p)
    else:
        return p


def main(args):
    # Setup logging
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

    # Create output directories if they don't exist
    logging.info("Check if output directory exists...")
    path = Path(args.output_dir)
    if not path.exists():
        logging.info("Output directory does not exist. Creating it...")
        path.mkdir(parents=True, exist_ok=True)
    else:
        logging.info("Output directory already exists.")

    logging.info(f"Current sample id: {args.sample_id}")
    logging.info("Load RNAseq data (10x format)...")
    if (os.path.exists(args.input_dir)):
        rnaseq_df = scprep.io.load_10X(parse_path(args.input_dir)).T
    else:
        raise Exception("Not a valid path for the input file")
    # Set genes as index
    rnaseq_df = rnaseq_df.rename_axis("index", axis=1)

    logging.info("Load database with CCIs...")
    lr_pairs = pd.read_csv(args.interactions_db)
    lr_pairs = lr_pairs.astype(str)

    # Metadata for the single cells
    logging.info("Load metadata...")
    meta_df = pd.read_csv(args.meta, index_col=0)
    meta_df.index.name = "index"
    meta_df = meta_df.loc[meta_df.index.isin(rnaseq_df.columns)]

    # Cell-cell Interactions and Communication Analysis
    # The pipeline integrates the RNA-seq and PPI datasets by using the analysis setups.
    # It generates an interaction space containing an instance for each sample/cell type, containing the values assigned to each protein in the PPI list given the setups for computing the CCI and CCC scores.

    logging.info("Infer CCIs")
    interactions = c2c.analysis.SingleCellInteractions(rnaseq_data=rnaseq_df,
                                                       ppi_data=lr_pairs,
                                                       metadata=meta_df,
                                                       interaction_columns=(
                                                           'source_genesymbol', 'target_genesymbol'),
                                                       communication_score='expression_gmean',
                                                       cci_score='bray_curtis',
                                                       cci_type='directed',
                                                       aggregation_method='average',
                                                       barcode_col='index',
                                                       celltype_col=args.annot,
                                                       complex_sep='_',
                                                       verbose=True)

    # **Compute communication scores for each PPI or LR pair**
    interactions.compute_pairwise_communication_scores()
    logging.info("Perform permutation analysis...")
    interactions.permute_cell_labels(evaluation='communication',
                                                permutations=args.n_perm,
                                                fdr_correction=True,
                                                verbose=True)
    logging.info("Save results...")
    interactions.ccc_permutation_pvalues.to_csv(
        f"{args.output_dir}/cell2cell__{args.sample_id}.csv", sep="\t")

    logging.info("Save interactions object as pickle...")
    with open(f"{args.output_dir}/cell2cell__{args.sample_id}.pickle", "wb") as handle:
        pickle.dump(interactions, handle, protocol=pickle.HIGHEST_PROTOCOL)

    logging.info("COMPLETED")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Infer CCIs with Cell2Cell")
    parser.add_argument("-f", "--input_dir", type=str,
                        help="Path to 10X data directory")
    parser.add_argument("-p", "--n_perm", type=int,
                        default=1000, help="Number of permutations for permutation testing")
    parser.add_argument("-db", "--interactions_db", type=str,
                        help="Path to custom database with interactions (csv)")
    parser.add_argument("-a", "--annot", type=str,
                        help="Column in metadata containing the cell type labels")
    parser.add_argument("-id", "--sample_id", type=str, default="")
    parser.add_argument("-o", "--output_dir", type=str,
                        default="output", help="Output directory")
    parser.add_argument("-m", "--meta", type=str,
                        help="Path to metadata file (CSV)")

    args = parser.parse_args()
    main(args)
