# Overview 

LIT PCBA docking is complete and being rsync'd to local from tobias (slow). 

The entries of intrest are the highest resolution PDB complexes within the set, located in: 

`./best_res_lit_pcba.csv` 

| Protein  | PDB_ID |
|----------|--------|
| ADRB2    | 4lde   |
| ALDH1    | 5l2m   |
| ESR1ago  | 2qzo   |
| ESR1ant  | 2iog   |
| FEN1     | 5fv7   |
| GBA      | 2v3d   |
| IDH1     | 4umx   |
| KAT2A    | 5mlj   |
| MAPK1    | 4zzn   |
| MTORC1   | 4dri   |
| OPRK1    | 6b73   |
| PKM2     | 3gr4   |
| PPARG    | 3b1m   |
| TP53     | 3zme   |
| VDR      | 3a2j   |

In simple csv: 

Protein,PDB_ID
ADRB2,4lde
ALDH1,5l2m
ESR1ago,2qzo
ESR1ant,2iog
FEN1,5fv7
GBA,2v3d
IDH1,4umx
KAT2A,5mlj
MAPK1,4zzn
MTORC1,4dri
OPRK1,6b73
PKM2,3gr4
PPARG,3b1m
TP53,3zme
VDR,3a2j

I will need a way to access the docking files related to these PDB_IDs/targets. 

# Docking Files 

Example with ADRB2: 

`ADRB2_4lde_active_glide_lib.sdfgz`
`ADRB2_4lde_inactive_glide_lib.sdfgz`

So, from the PDB ID of '4lde' for ADRB2, we would need to interact with the following files: 
- `ADRB2_4lde_active_glide_lib.sdfgz`
- `ADRB2_4lde_inactive_glide_lib.sdfgz`

In the directory: 

`/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2`

The steps we will need to do: 

1. Glide sort utility: keeping the best docking score per title. Unsure if we need to decompress from `.sdfgz` to `.sdf` first. However, we will need to decompress in order to run the strain calculation. 

2. Torsion Strain Calculation to produce the corresponding strain `.csv` file that will be operable with our papermill format. Papermill requires the path structures of the docking sdfs and the strain csv files, we will need to save the csv file in such a way as to be able to access it from the papermill notebook. The approach we used previously was: 

```py

import os
import pprint

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

# pprint.pprint(parameters_list)

# Execute the notebook for each set of parameters
for params in parameters_list:
    pm.execute_notebook(
        "gpcr_papermill.ipynb",
        params["output_notebook"],
        parameters=params["parameters"],
    )

```

First, we scan the directory for uppercase to define the subfolder we will process. The subfolder defines the title suffix is also used to find the desired files because the subfolder name is the same as the PDB ID (there was only one per target in this dataset). 

We will need a different approach.

We could still look in capital subfolder to define the `protein_name`, however that will need to be combined with our desired PDB IDs. Here is our directory structure: 

ADRB2/
ALDH1/
ESR1/
ESR1ago/
ESR1ant/
FEN1/
GBA/
IDH1/
KAT2A/
MAPK1/
MTORC1/
OPRK1/
PKM2/
PPARG/
TP53/
VDR/

However, we will need to build these file names: 

`ADRB2_4lde_active_glide_lib.sdfgz`
`ADRB2_4lde_inactive_glide_lib.sdfgz`

We could loop through a dictionary like structure like: 

```py
#key[value]
ADRB2[4lde]
```
The key will be used to search for the subfolder, so that would go to `ADRB2/`, then the value can be used to search for the file types 

(we will need to decompress the `.sdfgz` files to `.sdf` files first)

`key/key_value_active_glide_lib.sdf`
`key/key_value_inactive_glide_lib.sdf`

The strain value could be handled in a similar way, perhaps in a strain subfolder. Like this: 

`key/strain/key_value_active_glide_lib.csv`
`key/strain/key_value_inactive_glide_lib.csv`

We will need to build a dictionary like structure to handle this. We can probably do this from the protein pdb csv file directly. 

Let's continue to glide sort and torsion strain in the meanwhile, which will use a similar logic. 

# Glide Sort 

From history, I used this before: 

```sh
/opt/schrodinger/suites2023-3/utilities/glide_sort -o "/Users/lkv206/work/to_do_projects/chembl_ligands/GPCR-Bench-master/ADRB1/docking/ADRB1_active_docking_lib_sorted.sdf" -use_dscore -best_by_title "/Users/lkv206/work/to_do_projects/chembl_ligands/GPCR-Bench-master/ADRB1/docking/ADRB1_active_docking_lib.sdf"
```
Do I have an example with an sdfgz? No, not locally. The `glide_sort` utility mentions either sdfgz or sdf is fine for input, but states that it wants identical input/output file types, it may handle decompression however. 

Let's test: 

glide_sort -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted.sdf" -use_dscore -best_by_title "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_glide_lib.sdfgz"

Output: 

```
                        Glide Sort                                  
----------------------------------------------------------------------
REPORT OF BEST 33 POSES

The sorted ligand structures were written to the file:
    /Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted.sdf

Final rankings based on original docking score.

0 poses were rejected by the energy filters,
    Coul+vdw Energy    <=     0.0
    Hbond Interaction  <=     0.0
    Metal Interaction  <=    10.0
    GlideScore         <=   100.0
    Docking Score      <=   100.0
(If any of the above properties is not defined for a given pose,
 the corresponding filter is not applied to that pose.)
 ```

This seems good, though it is hard to tell without a before/after comparison. There also may not have been any duplicates. Let's also try the inactive file. 

```sh
glide_sort -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_inactive_docking_lib_sorted.sdf" -use_dscore -best_by_title "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_inactive_glide_lib.sdfgz"
```
```
Running glide_sort in large-file mode.

Glide Sort                                  
----------------------------------------------------------------------
REPORT OF BEST 103591 POSES

The sorted ligand structures were written to the file:
    /Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_inactive_docking_lib_sorted.sdf

Final rankings based on original docking score.

0 poses were rejected by the energy filters,
    Coul+vdw Energy    <=     0.0
    Hbond Interaction  <=     0.0
    Metal Interaction  <=    10.0
    GlideScore         <=   100.0
    Docking Score      <=   100.0
(If any of the above properties is not defined for a given pose,
 the corresponding filter is not applied to that pose.)

```
No poses rejected? 

This is a little odd. 

We should check for duplicates now, prior to the strain calculation. We will need an py/ipynb for this. 

Created `adrb2_duplicate_check.py` (interactive python style, psuedo-ipynb), it does appear to be to not have any duplicates. Odd, but possible. 

So both the conversion and deduplication seems to be fine here. We should commit here to update our progress and return to to strain. 