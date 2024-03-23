# %%
import pandas as pd
import os
from pprint import pprint
import subprocess

# print current working dir
print(os.getcwd())

# %%
best_res_lit_pcba = pd.read_csv("best_res_lit_pcba.csv")

# %%
best_res_lit_pcba

# %%
best_res_lit_pcba.columns

# %%
protein_pdbid = best_res_lit_pcba[["Protein", "PDB_ID"]]

# %%
# convert protein_pdbid dataframe to dictionary, where protein is the key and PDB_ID is the value
# we set the index of the dataframe to 'Protein', call 'PDB_ID' and then convert the pair to a dictionary
protein_pdbid_dict = protein_pdbid.set_index("Protein")["PDB_ID"].to_dict()
protein_pdbid_dict

# %%
file_id = []
for key, value in protein_pdbid_dict.items():
    file_id.append(f"{key}" + "-" + f"{value}")

print(file_id)

# %%
sdf_files = []

for dirpath, dirnames, filenames in os.walk("."):
    for filename in filenames:
        if any(id in filename for id in file_id) and filename.endswith(".sdf"):
            sdf_file_path = os.path.abspath(os.path.join(dirpath, filename))
            sdf_files.append(sdf_file_path)

print(sdf_files)
# %%
output_csv_files = []

for file in sdf_files:
    dir_name = os.path.dirname(file)  # Get the directory name
    base_name = os.path.basename(file)  # Get the base name
    new_dir_name = os.path.join(dir_name, 'strain')  # Add 'strain' to the directory
    csv_file = os.path.join(new_dir_name, base_name.replace(".sdf", ".csv"))  # Replace '.sdf' with '.csv' in the base name
    output_csv_files.append(csv_file)

print(output_csv_files)
# %%
print(sdf_files[0], output_csv_files[0])
# %%

# convert the sdfgz_files list and output_sdf_files list to a dictionary
sdf_csv_dict = dict(zip(sdf_files, output_csv_files))
print(sdf_csv_dict)
# %% [markdown]
# # Previously successful command:
# ```bash
# python ~/scripts/strain/refactor_Torsion_Strain.py -i "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_inactive_docking_lib_sorted.sdf" -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted_strain.csv"
# ```
# %%

# %%

for key, value in sdf_csv_dict.items():
    print(f"conda run analytics_env python ~/scripts/strain/refactor_Torsion_Strain.py -i {key} -o {value}")

# %%
commands = []
for key, value in sdf_csv_dict.items():
    commands.append(
        f"conda run analytics_env python ~/scripts/strain/refactor_Torsion_Strain.py -i {key} -o {value}"
    )
print(commands)
# %%

# to use parallel via GNU parallel, we need to split this list into new lines
commands_string = "\n".join(commands)
print(commands_string)
# %%
# set up the subprocess
process = subprocess.Popen(
    "parallel",
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
)
# %%

commands_string = "\n".join(commands)

# Setup subprocess to call parallel, sending the commands via stdin
process = subprocess.Popen(
    ["parallel"],  # Call GNU Parallel
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
)

# Pass the commands string to GNU Parallel through stdin
# and capture stdout and stderr
stdout, stderr = process.communicate(commands_string)

# handle the output and errors
print(stdout)
if stderr:
    print(f"Errors:\n{stderr}")

# saving stdout and stderr to files
with open("parallel_strain_writer.log", "w") as f_out, open(
    "parallel_errors.txt", "w"
) as f_err:
    f_out.write(stdout)
    if stderr:
        f_err.write(stderr)
