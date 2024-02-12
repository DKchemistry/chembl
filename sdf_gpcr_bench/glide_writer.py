import os
import sys

# Iterate over all subdirectories of the current directory
for dirpath, _, filenames in os.walk("."):
    depth = dirpath.count(os.sep)
    if depth == 1:
        os.makedirs(os.path.join(dirpath, 'docking'), exist_ok=True)
        # Get the zip and sdf files in the current subdirectory
        # need abs path for glide to work
        zip_files = [
            os.path.abspath(os.path.join(dirpath, file))
            for file in filenames
            if file.endswith(".zip")
        ]
        decoy_sdf = [
            os.path.abspath(os.path.join(dirpath, file))
            for file in filenames
            if file.endswith("decoy.sdf")
        ]
        active_sdf = [
            os.path.abspath(os.path.join(dirpath, file))
            for file in filenames
            if file.endswith("active.sdf")
        ]

        # If more than one criteria matches, exit the script
        if len(zip_files) > 1 or len(decoy_sdf) > 1 or len(active_sdf) > 1:
            print("More than one criteria matched. Exiting the script.")
            sys.exit()
        
        # if no criteria matches, print the directory and continue
        if len(zip_files) == 0 or len(decoy_sdf) == 0 or len(active_sdf) == 0:
            print(f"No criteria matched in {dirpath}.")
            continue

        # set names to be used in file writing
        protein_name = os.path.basename(dirpath)  # ex ADRB1

        decoy_docking_in = os.path.join(
            dirpath, "docking", f"{protein_name}_decoy_docking.in"
        )
        # print(decoy_docking_in)

        active_docking_in = os.path.join(
            dirpath, "docking", f"{protein_name}_active_docking.in"
        )
        # print(active_docking_in)

        decoy_docking_sh = os.path.join(
            dirpath, "docking", f"{protein_name}_decoy_docking.sh"
        )
        # print(decoy_docking_sh)

        active_docking_sh = os.path.join(
            dirpath, "docking", f"{protein_name}_active_docking.sh"
        )
        # print(active_docking_sh)

        with open(decoy_docking_in, "w") as f:
            f.write(f"GRIDFILE {zip_files[0]}\n")
            f.write(f"LIGANDFILE {decoy_sdf[0]}\n")
            f.write("POSE_OUTTYPE ligandlib_sd\n")
            f.write("DOCKING_METHOD confgen\n")
            f.write("PRECISION SP\n")
            f.write("AMIDE_MODE penal\n")
            f.write("SAMPLE_RINGS True\n")
            f.write("EPIK_PENALTIES True\n")

        with open(active_docking_in, "w") as f:
            f.write(f"GRIDFILE {zip_files[0]}\n")
            f.write(f"LIGANDFILE {active_sdf[0]}\n")
            f.write("POSE_OUTTYPE ligandlib_sd\n")
            f.write("DOCKING_METHOD confgen\n")
            f.write("PRECISION SP\n")
            f.write("AMIDE_MODE penal\n")
            f.write("SAMPLE_RINGS True\n")
            f.write("EPIK_PENALTIES True\n")

        with open(decoy_docking_sh, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write(
                f"$SCHRODINGER/glide -HOST localhost:10 -NJOBS 20 -OVERWRITE -JOBNAME {protein_name}_decoy_docking {protein_name}_decoy_docking.in"
            )
        os.chmod(decoy_docking_sh, 0o755)

        with open(active_docking_sh, "w") as f:
            f.write("#!/bin/bash\n\n")
            f.write(
                f"$SCHRODINGER/glide -HOST localhost:2 -NJOBS 4 -OVERWRITE -JOBNAME {protein_name}_active_docking {protein_name}_active_docking.in"
            )
        os.chmod(active_docking_sh, 0o755)
