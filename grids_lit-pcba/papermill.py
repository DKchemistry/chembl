# %%
import papermill as pm
import os
import pprint
import pandas as pd

# %%
# vs code is so fucking dumb
# have to remove the following code later
os.chdir("/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba")

print(os.getcwd())

# %%
pm.inspect_notebook("litpcba_papermill.ipynb")


# %%
# # execute_notebook(<input notebook>, <output notebook>, <dictionary of parameters>)

# pm.execute_notebook(
#     "gpcr_papermill.ipynb",
#     "./papermill/notebooks/gpcr_papermill_output.ipynb",
#     parameters={
#         "title_suffix": "GPR40",
#         "file_path_sdf_active": "./GPR40/docking/GPR40_active_docking_lib_sorted.sdf",
#         "file_path_sdf_decoy": "./GPR40/docking/GPR40_decoy_docking_lib_sorted.sdf",
#         "file_path_strain_active": "./GPR40/strain/GPR40_active_docking_lib_sorted.csv",
#         "file_path_strain_decoy": "./GPR40/strain/GPR40_decoy_docking_lib_sorted.csv",
#     },
# )


# %%
# parameters = {
#     "title_suffix": "GPR40",
#     "file_path_sdf_active": "./GPR40/docking/GPR40_active_docking_lib_sorted.sdf",
#     "file_path_sdf_decoy": "./GPR40/docking/GPR40_decoy_docking_lib_sorted.sdf",
#     "file_path_strain_active": "./GPR40/strain/GPR40_active_docking_lib_sorted.csv",
#     "file_path_strain_decoy": "./GPR40/strain/GPR40_decoy_docking_lib_sorted.csv",
# }

# output_notebook = (
#     f"./papermill/notebooks/gpcr_papermill_output_{parameters['title_suffix']}.ipynb"
# )

# pm.execute_notebook(
#     "gpcr_papermill.ipynb",
#     output_notebook,
#     parameters=parameters,
# )

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

parameters_list = []
parameters_dict = {}

for id in file_id:
    parameters_list = []
    parameters_list.append(
        os.path.abspath(
            os.path.join("papermill", "notebooks", (id + "_litpcba_papermill.ipynb"))
        )
    )
    # print(parameters_list)

    for dirpath, dirname, filenames in os.walk("."):
        for filename in filenames:
            if filename.endswith(f"{id}" + "_active_glide_lib_sorted.sdf"):
                parameters_list.append(os.path.abspath(filename))
            if filename.endswith(f"{id}" + "_inactive_glide_lib_sorted.sdf"):
                parameters_list.append(os.path.abspath(filename))
            if filename.endswith(f"{id}" + "_active_glide_lib_sorted.csv"):
                parameters_list.append(
                    os.path.abspath(os.path.join("strain", filename))
                )
            if filename.endswith(f"{id}" + "_inactive_glide_lib_sorted.csv"):
                parameters_list.append(
                    os.path.abspath(os.path.join("strain", filename))
                )
            parameters_dict[id] = parameters_list

for key, value in parameters_dict.items():
    print(key)
    for items in value:
        print(items)
    print("\n")

# %%

parameters_list = []
parameters_dict = {}

for id in file_id:
    # set the key to the id and the value to an empty dictionary
    parameters_dict[id] = {}
    parameters_dict[id]["title_suffix"] = id
    parameters_dict[id]["output_notebook"] = os.path.abspath(
        os.path.join("papermill", "notebooks", (id + "_litpcba_papermill.ipynb"))
    )

    for dirpath, dirname, filenames in os.walk("."):
        for filename in filenames:
            if filename.endswith(f"{id}" + "_active_glide_lib_sorted.sdf"):
                parameters_dict[id]["file_path_sdf_active"]=(os.path.abspath(filename))
            if filename.endswith(f"{id}" + "_inactive_glide_lib_sorted.sdf"):
                parameters_dict[id]["file_path_sdf_inactive"] = os.path.abspath(filename)
            if filename.endswith(f"{id}" + "_active_glide_lib_sorted.csv"):
                parameters_dict[id]["file_path_strain_active"] = os.path.abspath(
                    os.path.join("strain", filename)
                )
            if filename.endswith(f"{id}" + "_inactive_glide_lib_sorted.csv"):
                parameters_dict[id]["file_path_strain_inactive:"] = os.path.abspath(
                    os.path.join("strain", filename)
                )

for key, value in parameters_dict.items():
    print(key)
    for key, value in value.items():
        print(key, value)
    print("\n")
# %%

parameters_list = []

for id in file_id:
    parameters_dict = {}

    parameters_dict["output_notebook"] = os.path.abspath(
        os.path.join("papermill", "notebooks", (id + "_litpcba_papermill.ipynb"))
    )
    # my mistake was not initalizing the empty dictionary first
    parameters_dict["parameters"] = {}

    parameters_dict["parameters"]["title_suffix"] = id

    for dirpath, dirname, filenames in os.walk("."):
        for filename in filenames:
            if filename.endswith(f"{id}" + "_active_glide_lib_sorted.sdf"):
                parameters_dict["parameters"]["file_path_sdf_active"] = os.path.abspath(
                    filename
                )
            if filename.endswith(f"{id}" + "_inactive_glide_lib_sorted.sdf"):
                parameters_dict["parameters"]["file_path_sdf_decoy"] = (
                    os.path.abspath(filename)
                )
            if filename.endswith(f"{id}" + "_active_glide_lib_sorted.csv"):
                parameters_dict["parameters"]["file_path_strain_active"] = (
                    os.path.abspath(os.path.join("strain", filename))
                )
            if filename.endswith(f"{id}" + "_inactive_glide_lib_sorted.csv"):
                parameters_dict["parameters"]["file_path_strain_decoy"] = (
                    os.path.abspath(os.path.join("strain", filename))
                )
    # print(parameters_dict)
    # print("\n")
    parameters_list.append(parameters_dict)
pprint.pprint(parameters_list)
print(len(parameters_list))
# print(parameters_list[0])
# %%

for parameter_key, parameter_value in parameters_dict.items():
    print(parameter_key, parameter_value)
    print("litpcba_papermill.ipynb", parameter_key["output_notebook"], parameter_value)
    # pm.execute_notebook("litpcba_papermill.ipynb", parameter_value["output_notebook"], parameter_value)

print(parameters_dict[id])
for value in parameters_dict[id].values():
    print(value)
# %%

# %%

# Get a list of all subfolders in the current working directory that start with a capital letter
subfolders = [f.name for f in os.scandir(".") if f.is_dir() and f.name[0].isupper()]

parameters_list = []

# Create a parameters dictionary for each subfolder
for subfolder in subfolders:
    parameters = {
        "title_suffix": subfolder,
        "file_path_sdf_active": f"./{subfolder}/docking/{subfolder}_active_docking_lib_sorted.sdf",
        "file_path_sdf_decoy": f"./{subfolder}/docking/{subfolder}_decoy_docking_lib_sorted.sdf",
        "file_path_strain_active": f"./{subfolder}/strain/{subfolder}_active_docking_lib_sorted.csv",
        "file_path_strain_decoy": f"./{subfolder}/strain/{subfolder}_decoy_docking_lib_sorted.csv",
    }

    output_notebook = f"./papermill/notebooks/gpcr_papermill_output_{parameters['title_suffix']}.ipynb"

    parameters_list.append(
        {
            "output_notebook": output_notebook,
            "parameters": parameters,
        }
    )


pprint.pprint(parameters_list)

# # Execute the notebook for each set of parameters
# for params in parameters_list:
#     pm.execute_notebook(
#         "gpcr_papermill.ipynb",
#         params["output_notebook"],
#         parameters=params["parameters"],
#     )

# %%
pprint.pprint(parameters_list)

# %%
