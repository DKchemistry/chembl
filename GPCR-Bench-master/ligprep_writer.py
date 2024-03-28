import os
import glob
import argparse


def find_smi_files():
    """
    Search for .smi files within '*/stereochemistry/*'.
    """
    return glob.glob("*/stereochemistry/*.smi")


def write_inp_file(directory):
    """
    Write or preview the ligprep_noepik.inp file with fixed content.
    """
    content = """MAX_ATOMS   500
FORCE_FIELD   16
IONIZATION   0
USE_DESALTER   no
GENERATE_TAUTOMERS   no
DETERMINE_CHIRALITIES   no
IGNORE_CHIRALITIES   no
NUM_STEREOISOMERS   1
"""
    inp_file_path = os.path.join(directory, "ligprep_noepik.inp")
    if args.run:
        with open(inp_file_path, "w") as file:
            file.write(content)
    else:
        print(f"Preview for {inp_file_path}:\n{content}")


def write_sh_file(smi_file):
    """
    Generate and preview/write the .sh file for the given .smi file.
    """
    directory = os.path.dirname(smi_file)
    smi_filename = os.path.basename(smi_file)
    base_name = smi_filename.replace("_sc.smi", "")
    active = "active" in smi_filename
    njobs = "1" if active else "5"
    ncpus = "1" if active else "5"
    jobname = base_name + "_ligprep"
    sh_content = f"""#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS {njobs} -JOBNAME {jobname} -HOST localhost:{ncpus} -ismi {smi_filename} -osd {base_name}.sdf
"""
    sh_file_path = os.path.join(directory, base_name + "_ligprep.sh")
    if args.run:
        with open(sh_file_path, "w") as file:
            file.write(sh_content)
        os.chmod(sh_file_path, 0o755)
    else:
        print(f"Preview for {sh_file_path}:\n{sh_content}")


def main():
    smi_files = find_smi_files()
    for smi_file in smi_files:
        directory = os.path.dirname(smi_file)
        write_inp_file(directory)
        write_sh_file(smi_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automate the generation of .inp and .sh files for stereochemistry analysis."
    )
    parser.add_argument(
        "-run",
        action="store_true",
        help="Actually write the files. Without this, the script will run in dry mode.",
    )
    args = parser.parse_args()

    main()
