import argparse
import logging
from cellphonedb.utils import db_utils
import os
import fnmatch


def main(args):
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
    db_utils.create_db(args.input_dir)

    # filename = fnmatch.filter(os.listdir(args.input_dir), '*.zip')[-1]
    # os.rename(f"{args.input_dir}/{filename}", f"{args.input_dir}/cpdb.zip")
    logging.info("COMPLETED")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create CellPhoneDB database from files")

    # Parse arguments
    # Always required
    parser.add_argument("-i", "--input_dir", type=str,
                        default="", help="Directory wit files necessary to create the database")
    args = parser.parse_args()
    # TODO make sure that the line is commented when running from the command line, only uncomment if run from Python
    # args.input_dir = "/Users/joankant/Desktop/gaitigroup/Users/Joan/GBM_CCI_Analysis/output/wip_interactions_db"
    main(args)
