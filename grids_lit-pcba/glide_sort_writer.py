# %%
import pandas as pd
import os
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
protein_pdbid = best_res_lit_pcba[['Protein', 'PDB_ID']]

# %%
# convert protein_pdbid dataframe to dictionary, where protein is the key and PDB_ID is the value
# we set the index of the dataframe to 'Protein', call 'PDB_ID' and then convert the pair to a dictionary
protein_pdbid_dict = protein_pdbid.set_index('Protein')['PDB_ID'].to_dict()
protein_pdbid_dict

# %%
file_id = []
for key, value in protein_pdbid_dict.items():
    file_id.append(f"{key}" + "_" + f"{value}")

print(file_id)

# %%
sdfgz_files = []

for dirpath, dirnames, filenames in os.walk("."):
  for filename in filenames:
    if any(id in filename for id in file_id) and filename.endswith(".sdfgz"):
      sdfgz_file_path = os.path.abspath(os.path.join(dirpath, filename))
      sdfgz_files.append(sdfgz_file_path)

print(sdfgz_files)
# %%
output_sdf_files = []

for file in sdfgz_files:
  sdf_file = file.replace(".sdfgz", ".sdf")
  output_sdf_files.append(sdf_file)

print(output_sdf_files)
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

for key, value in sdfgz_sdf_dict.items():
  print(f"glide_sort -o {value} -use_dscore -best_by_title {key}")

# %%
