# %%
import papermill as pm
import os
import pprint
import pandas as pd

from ploomber import DAG
from ploomber.products import File
from ploomber.tasks import NotebookRunner
from ploomber.executors import Parallel
from pathlib import Path

# %%

SOURCE_NOTEBOOK = os.path.abspath("litpcba_papermill.ipynb")
print(SOURCE_NOTEBOOK)

# %%
print(os.getcwd())


# %%
pm.inspect_notebook(SOURCE_NOTEBOOK)

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
# Get the current working directory (base path)
base_path = os.getcwd()

parameters_list = []

for id in file_id:
    protein_name, _ = id.split("-")  # Split the id to extract the protein name.

    parameters_dict = {}
    # Construct the path for output_notebook relative to base_path
    parameters_dict["output_notebook"] = os.path.join(
        base_path, "ploomer_output", "notebooks", (id + "_litpcba_ploomer.ipynb")
    )
    parameters_dict["parameters"] = {}
    parameters_dict["parameters"]["title_suffix"] = id

    for dirpath, dirnames, filenames in os.walk(base_path):
        for filename in filenames:
            # Construct each path by checking if the filename ends with the required suffix
            if filename.endswith(f"{id}_active_glide_lib_sorted.sdf"):
                # Join the path without repeating the base path
                parameters_dict["parameters"]["file_path_sdf_active"] = os.path.join(
                    dirpath, filename
                )
            if filename.endswith(f"{id}_inactive_glide_lib_sorted.sdf"):
                parameters_dict["parameters"]["file_path_sdf_decoy"] = os.path.join(
                    dirpath, filename
                )
            if filename.endswith(f"{id}_active_glide_lib_sorted.csv"):
                parameters_dict["parameters"]["file_path_strain_active"] = os.path.join(
                    dirpath, filename
                )
            if filename.endswith(f"{id}_inactive_glide_lib_sorted.csv"):
                parameters_dict["parameters"]["file_path_strain_decoy"] = os.path.join(
                    dirpath, filename
                )
    parameters_list.append(parameters_dict)

# pprint.pprint(parameters_list)
pprint.pprint(parameters_list)
print(len(parameters_list))
# print(parameters_list[0])
# %%

for params in parameters_list:
    print(SOURCE_NOTEBOOK, params["output_notebook"], params["parameters"])

# %%
# Define your DAG with a Parallel executor for parallel execution
dag = DAG(executor=Parallel())

for params in parameters_list:
    output_notebook_path = params["output_notebook"]

    # Setup a task for running the notebook
    # Note: No need to specify 'papermill_params' or 'engine_name'
    NotebookRunner(
        Path(SOURCE_NOTEBOOK),
        File(output_notebook_path),
        dag=dag,
        params=params["parameters"],
    )

if __name__ == "__main__":
    dag.build(force=True)

# %%
