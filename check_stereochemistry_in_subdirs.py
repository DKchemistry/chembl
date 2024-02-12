import os
import argparse
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import rdchem


def check_stereochemistry(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    has_unspecified_stereo = any(not assigned for center, assigned in chiral_centers)

    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            if bond.GetStereo() <= rdchem.BondStereo.STEREOANY:
                has_unspecified_stereo = True
                break

    return not has_unspecified_stereo, ""


def process_smi_files(directory):
    log = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".smi"):
                file_path = os.path.join(root, file)
                with open(file_path, "r") as f:
                    for line in f:
                        if (
                            line.lower()
                            .strip()
                            .startswith(
                                ("smiles", '"smiles"', "smiles id", '"smiles" "id"')
                            )
                        ):
                            continue
                        parts = line.strip().split()
                        if len(parts) < 2:  # Ensure there are at least two parts
                            continue
                        smiles, id_ = parts[:2]
                        is_fully_enumerated, reason = check_stereochemistry(smiles)
                        if not is_fully_enumerated:
                            log_entry = f"{file_path}: {smiles} with ID {id_} is not fully enumerated. Reason: {reason}"
                            print(log_entry)
                            log.append(log_entry)

    # Get the current time for timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"enumeration_log_{timestamp}.txt"
    with open(log_filename, "w") as log_file:
        for entry in log:
            log_file.write(entry + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Check if SMILES strings in .smi files within subdirectories have fully enumerated stereochemistry. Example usage: python script.py -d '/path/to/directory'"
    )
    parser.add_argument(
        "-d",
        "--top-directory",
        dest="top_directory",
        type=str,
        required=True,
        help="Short flag for the top-level directory to search for .smi files. Example: '/path/to/directory'",
    )

    args = parser.parse_args()
    process_smi_files(args.top_directory)


if __name__ == "__main__":
    main()
