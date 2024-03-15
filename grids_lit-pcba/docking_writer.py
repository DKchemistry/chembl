import os
import glob

# Define the base directory where to look for subdirectories
base_dir = "/mnt/data/dk/work/lit-pcba"

# Walk through the current directory
for root, dirs, files in os.walk("."):
    for file in files:
        # Skip files that do not end with "_grid.zip"
        if not file.endswith("_grid.zip"):
            continue

        prefix = file.split("_")[0]  # Get prefix up to the first '_'
        name = "_".join(file.split("_")[:2])  # Get name up to the second '_'

        # Match subdirectories in the base directory
        matched_dir = None
        for dir in glob.glob(os.path.join(base_dir, prefix + "*")):
            if os.path.isdir(dir):
                matched_dir = dir
                break

        if matched_dir:
            # Construct file paths for actives and inactives
            actives_ligand_abs_path = os.path.join(
                matched_dir, "actives_rdkit_ligprep.sdf"
            )
            inactives_ligand_abs_path = os.path.join(matched_dir, "inactives_rdkit.sdf")

            # Check if the ligand files exist
            if os.path.exists(actives_ligand_abs_path) and os.path.exists(
                inactives_ligand_abs_path
            ):
                print(f"Found matching directory: {matched_dir}")
                # Write input files
                for ligand_type in ["active", "inactive"]:
                    in_file_path = os.path.join(root, f"{name}_{ligand_type}.in")
                    ligand_abs_path = (
                        actives_ligand_abs_path
                        if ligand_type == "active"
                        else inactives_ligand_abs_path
                    )
                    with open(in_file_path, "w") as in_file:
                        in_file.write(
                            f"GRIDFILE {os.path.abspath(os.path.join(root, file))}\n"
                            f"LIGANDFILE {ligand_abs_path}\n"
                            "POSE_OUTTYPE ligandlib_sd\n"
                            "DOCKING_METHOD confgen\n"
                            "PRECISION SP\n"
                            "AMIDE_MODE penal\n"
                            "SAMPLE_RINGS True\n"
                            "EPIK_PENALTIES True\n"
                        )
                    print(f"Written {in_file_path}")

                # Write shell scripts
                for ligand_type in ["active", "inactive"]:
                    sh_file_path = os.path.join(root, f"{name}_{ligand_type}.sh")
                    with open(sh_file_path, "w") as sh_file:
                        njobs = "1" if ligand_type == "active" else "450"
                        host = (
                            "localhost:1"
                            if ligand_type == "active"
                            else "localhost:100"
                        )
                        sh_file.write(
                            f"/mnt/data/dk/Schrodinger_adv_2021_1/glide -HOST {host} -NJOBS {njobs} "
                            f"-OVERWRITE -JOBNAME {name}_{ligand_type}_glide {name}_{ligand_type}.in\n"
                        )
                    print(f"Written {sh_file_path}")
            else:
                print(f"Ligand files not found in {matched_dir}")
        else:
            print(f"No matching directory found for {file}")

print("Script execution completed.")
