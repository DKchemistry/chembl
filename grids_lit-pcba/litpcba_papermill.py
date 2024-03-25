# %%
import papermill as pm
import os
import pprint
import pandas as pd

# %%
print(os.getcwd())

# %%
pm.inspect_notebook("litpcba_papermill.ipynb")

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
        base_path, "papermill", "notebooks", (id + "_litpcba_papermill.ipynb")
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
    print("litpcba_papermill.ipynb", params["output_notebook"], params["parameters"])
    # pm.execute_notebook("litpcba_papermill.ipynb", parameter_value["output_notebook"], parameter_value)

# %%

# %%

# # Get a list of all subfolders in the current working directory that start with a capital letter
# subfolders = [f.name for f in os.scandir(".") if f.is_dir() and f.name[0].isupper()]

# parameters_list = []

# # Create a parameters dictionary for each subfolder
# for subfolder in subfolders:
#     parameters = {
#         "title_suffix": subfolder,
#         "file_path_sdf_active": f"./{subfolder}/docking/{subfolder}_active_docking_lib_sorted.sdf",
#         "file_path_sdf_decoy": f"./{subfolder}/docking/{subfolder}_decoy_docking_lib_sorted.sdf",
#         "file_path_strain_active": f"./{subfolder}/strain/{subfolder}_active_docking_lib_sorted.csv",
#         "file_path_strain_decoy": f"./{subfolder}/strain/{subfolder}_decoy_docking_lib_sorted.csv",
#     }

#     output_notebook = f"./papermill/notebooks/gpcr_papermill_output_{parameters['title_suffix']}.ipynb"

#     parameters_list.append(
#         {
#             "output_notebook": output_notebook,
#             "parameters": parameters,
#         }
#     )


# pprint.pprint(parameters_list)

# # Execute the notebook for each set of parameters
# for params in parameters_list:
#     pm.execute_notebook(
#         "gpcr_papermill.ipynb",
#         params["output_notebook"],
#         parameters=params["parameters"],
#     )
# %%
