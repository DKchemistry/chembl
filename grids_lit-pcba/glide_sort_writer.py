# %%
import pandas as pd
import os
import subprocess
from pprint import pprint


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
    file_id.append(f"{key}" + "_" + f"{value}")

print(file_id)

# %%
sdfgz_files = []
output_sdf_files = []

for dirpath, dirnames, filenames in os.walk("."):
    for filename in filenames:
        if any(id in filename for id in file_id) and filename.endswith(".sdfgz"):
            sdfgz_file_path = os.path.abspath(os.path.join(dirpath, filename))
            print("Found sdfgz input file:")
            print(sdfgz_file_path)
            sdfgz_files.append(sdfgz_file_path)
            # define output path naming convention
            parts = os.path.basename(sdfgz_file_path).split("_")
            parts_join = "-".join(parts[:2])+"_"+"_".join(parts[2:])
            sdf_file = os.path.abspath(os.path.join(
                dirpath, parts_join.replace(".sdfgz", "_sorted.sdf")
            ))
            print("Found sorted sdf output file:")
            print(sdf_file)
            output_sdf_files.append(sdf_file)

# %%
print(sdfgz_files[0], output_sdf_files[0])
# %%

# convert the sdfgz_files list and output_sdf_files list to a dictionary
sdfgz_sdf_dict = dict(zip(sdfgz_files, output_sdf_files))
print(sdfgz_sdf_dict)
# %% [markdown]
# # Previously successful command:
# ```sh
# glide_sort -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted.sdf" -use_dscore -best_by_title "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_glide_lib.sdfgz"
# ```
# %%
#import pyperclip

# %%

for key, value in sdfgz_sdf_dict.items():
    print(f"glide_sort -o {value} -use_dscore -best_by_title {key}")

# %%

# copy = []
# for key, value in sdfgz_sdf_dict.items():
# #  print(f"glide_sort -o {value} -use_dscore -best_by_title {key}")
#   copy.append(f"glide_sort -o {value} -use_dscore -best_by_title {key}")

# pyperclip.copy("\n".join(copy))
# pyperclip.paste()
# %%
commands = []
for key, value in sdfgz_sdf_dict.items():
    commands.append(f"glide_sort -o {value} -use_dscore -best_by_title {key}")
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
# pass the commands to the subprocess stdin via communicate, the passing to stdin is implicit, you can't set communicate(stdin=commands_string).
stdout, stderr = process.communicate(commands_string)

# print the output
print(stdout)
if stderr:
    print(f"Errors:\n{stderr}")

# Save stdout and stderr to a file
with open('glide_sort_writer_output.txt', 'w') as f:
    f.write("STDOUT:\n")
    f.write(stdout)
    if stderr:
        f.write("\nSTDERR:\n")
        f.write(stderr)
# %%
