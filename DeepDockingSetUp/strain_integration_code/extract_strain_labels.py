import argparse
import builtins as __builtin__
import glob
import gzip
import os
from contextlib import closing
from multiprocessing import Pool


# For debugging purposes only:
# def print(*args, **kwargs):
#     __builtin__.print("\t extract_strain: ", end="")
#     return __builtin__.print(*args, **kwargs)

TEST_PATH = "/Users/lkv206/work/to_do_projects/chembl_ligands/GPCR-Bench-master/5HT1B/strain/5HT1B_active_docking_lib_sorted.csv"
TEST_OUTPUT = "/Users/lkv206/work/to_do_projects/chembl_ligands/DeepDockingSetUp/strain_labels/extract_strain_labels.csv"

parser = argparse.ArgumentParser(
    description="Extract strain labels from a set of images"
)
parser.add_argument(
    "-cd",
    "--csv_directory",
    type=str,
    help="Directory containing the csv files",
    required=False,
    default=TEST_PATH,
)
parser.add_argument(
    "-o",
    "--output_filepath",
    type=str,
    help="Output file",
    required=False,
    default="/Users/lkv206/work/to_do_projects/chembl_ligands/DeepDockingSetUp/strain_labels/extract_strain_labels.csv",
)

# When coding in an interactive env
# It will SystemExit:2, so we need to catch it
try:
    args = parser.parse_args()
    csv_directory = args.csv_directory
    output_filepath = args.output_filepath
except SystemExit:
    csv_directory = TEST_PATH
    output_filepath = TEST_OUTPUT

# # Comment these lines back in when done testing
# args = parser.parse_args()
# csv_directory = args.csv_directory

print(f"csv_directory: {csv_directory}")
print(f"output_filepath: {output_filepath}")


def extract_strain_labels(csv_file):
    strain_labels = []
    with open(csv_file, "r") as f:
        for line in f:
            id = str(line.split(",")[0])
            strain = float(line.split(",")[1])
            if strain < 0:
                continue
            if strain > 0:
                # print(f"Name: {id}, Strain: {strain}")
                strain_labels.append([strain, id])
    return strain_labels


if __name__ == "__main__":
    results = extract_strain_labels(csv_directory)
    print(results)

    with open(output_filepath, "w") as f:
        for result in results:
            f.write(f"{result[0]},{result[1]}\n")

    print(f"Strain labels saved to {output_filepath}")
