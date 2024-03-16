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

```sh
glide_sort -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted.sdf" -use_dscore -best_by_title "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_glide_lib.sdfgz"
```

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

Commit has been added. Rsync is still in progress (quite slow right now). 

When running strain on Tobias, I had to do some form of refactoring in order to successfully import the XML library, I do not believe that function is local to this mac. Let's go look on tobias for what I did. 

So, the change seems work like so: 

In `refactor_TL_Functions.py`

```py 

import xml.etree.ElementTree as ET
from rdkit import Chem
import os
import numpy as np
from math import sqrt, atan2, pi
from math import ceil

#tree = ET.parse("TL_2.1_VERSION_6.xml")
#root = tree.getroot()
# Determine the directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))

# Construct the absolute path to the XML file
xml_file_path = os.path.join(script_dir, "TL_2.1_VERSION_6.xml")

# Use the absolute path to parse the XML file
tree = ET.parse(xml_file_path)
root = tree.getroot()

```
This is contained in a directory in `/mnt/data/dk/scripts/strain` as such: 

```sh 
TL_2.1_VERSION_6.xml
__pycache__
refactor_TL_Functions.py
refactor_Torsion_Strain.py
```

I should rsync to a local strain calculation. I should also check my history of how this command was run and where it is on my path on tobias. 

The path (tobias) seems to be as expected: 

`export PATH=/mnt/data/dk/scripts:$PATH`

The command was run via a `strain_runner.sh` script that required some args and did a dry run prior to execution. 

```sh 
./strain_runner.sh -b . -s /mnt/data/dk/scripts/strain/refactor_Torsion_Strain.py -e rdkit_en
```
-b is for the directory containing the subdirectories, -s is the path to the strain (as we can see), and -e was the conda environment. It does a confirmation step so we can see what the commands would be prior to execution (I used parallel here). 

```sh
python /mnt/data/dk/scripts/strain/refactor_Torsion_Strain.py -i "/mnt/data/dk/work/gpcr_bench/GPCR-Bench-master/SMO/docking/SMO_decoy_docking_lib.sdf" -o "/mnt/data/dk/work/gpcr_bench/GPCR-Bench-master/SMO/strain/SMO_decoy_docking_lib.csv"
```

The strain runner (and the sort runner), will probably need to be refactored to handle the new directory structure. 

First, let's get the local version of our script running because the latency to tobias is pretty rough. It appears I already have it. So we can start testing the strain. 

# Strain Calculation 

One caveat with our python strain script is that despite exporting it to path, it can not be easily called from the command line as we need to pass 'python' prior. A shebang can be added to the top of the script to make it executable, but I am not sure if that is the best approach as we would need the conda environment with rdkit - reducing portability of the script and adding later confusion. It is best to call the full path of the script for now. 

```sh
python ~/scripts/strain/refactor_Torsion_Strain.py -i "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted.sdf" -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/strain/ADRB2_4lde_active_docking_lib_sorted_strain.csv"
```

Unfortunately, the script can not natively force the creation of the strain directory I specified in the output. 

I refactored and committed the refactor to allow for creating the dir. 

Output: 

```
33 molecules finished reading. Calculating strain energy...
33 successful / 0 NA. Please check: /Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/strain/ADRB2_4lde_active_docking_lib_sorted_strain.csv

```
Now let's run the inactive if it computes and continue testing the duplicates. 

```sh
python ~/scripts/strain/refactor_Torsion_Strain.py -i "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_inactive_docking_lib_sorted.sdf" -o "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted_strain.csv"
```

Seems to be running fine, it's just gonna take awhile. Once we have tested it for duplicates, we can automate it and run the rest of the targets, local will be fine. 

Output: 

```
103591 molecules finished reading. Calculating strain energy...
102674 successful / 917 NA. Please check: /Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/ADRB2/ADRB2_4lde_active_docking_lib_sorted_strain.csv

```

NOTE: Mistake in output file path, saved it as "active" and did not save it to strain. Manually `mv` to fix it since no overwrite occured. 

# Strain Metrics 

Initial tests using the adrb2 system look okay, but some things should be fixed and potentially retrospectively fixed for the GPCR-bench dataset. 

1. Even though we should only have `102674 successful / 917 NA.` molecules in the inactives, it appears as though we have 103K. This needs to be explored by examining the difference (if any) between the input csv of strain energy, and what the dataframe actually looks like. 

* Follow Up: 

```py
decoy_strain['Total_E'].isna().sum()
``` 
Returns 917. These are NaN values. This makes sense, however they will need to be dropped.

We can check the actives as well: 

```py
active_strain['Total_E'].isna().sum()
```
Returns 0. We should still drop these in case other active strain files have NaN values. 

2. Some strain values are flagged (negative energies) and should be removed. 

* Follow Up: 

```py
(decoy_strain['Total_E'] < 0).sum()
```
Returns 146. These are negative values. They will need to be dropped. 

We can check the actives as well: 

```py
(active_strain['Total_E'] < 0).sum()
```
Returns 0. We should still drop these in case other active strain files have negative values. 

Still, this doesn't make sense. I can see from `plot_density_strain` that there are negative values. Something is wrong. Is this a consequence of the density plot? It appears to be so. 

See here, if we run: 

```py 
plot_densities_strain(all_data, title_suffix)
```
We get a density plot with negative values for the actives: 

![alt text](image-1.png)

Yet, if we run: 

```py
all_data[(all_data['Total_E']>0) & (all_data['Activity']==1)]['Total_E'].hist()
```

We get a histogram with no negative values: 

![alt text](image-2.png)

By changing the function to make histograms instead (but normalized to an area of 1) with bins=50, like so below, we can see more clearly what is going on. The inactives do have some negative values, but the actives do not. 

```py
def plot_histogram_strain(df, title_suffix):
    # Hardcoded column names
    activity_col = "Activity"
    score_col = "Total_E"
    
    plt.hist(df.loc[df[activity_col] == 0, score_col], bins=50, label="Inactive", alpha=0.5, density=True)
    plt.hist(df.loc[df[activity_col] == 1, score_col], bins=50, label="Active", alpha=0.5, density=True)

    
    # Add title and labels
    plt.title(f"Histogram of Strain Energy for Active and Decoy Molecules ({title_suffix})")
    plt.xlabel("Total Strain Energy")
    plt.ylabel("Density")
    plt.legend(loc="best")
    
    # Show the plot
    plt.show()
```
```py
plot_histogram_strain(all_data, title_suffix)
```

![alt text](image-3.png)

This was perhaps a long detour to better understand how KDE plots can work and occasionally be deceptive compared to histograms. 

In conclusion, we still should still drop the NaN and negative values from the decoy data prior to the merging or concat steps.

We will need to integrate these changes into GPCR-bench analysis workflow as well and update the metrics. 

3. There appears to be some sort of issue in reading in the csv file containing the strain energies. For the ADRB2 system, I am seeing 33/33 successful strain calculations yet `active_strain` has 32 rows. 

* Follow Up: 

(brief note, remember our `cmc` function, `cmc <some command>` will copy the command and output)

```sh
$ wc -l ADRB2_4lde_active_docking_lib_sorted_strain.csv
      33 ADRB2_4lde_active_docking_lib_sorted_strain.csv
```
So the CSV really does have 33 entries. Why is the first one being lost?

The issue was with this function: 

```py
concatenate_csv_files(file_list)
```
I used header=0, where I should have used header=None. I believe header=0 specified that the first row was the headers, even though the function does rename them with the `column_names` variable. This will also need to be updated in the GPCR-bench workflow. 

I have added the commits to the repo for this file and `adrb2_duplicate_check.py`. 

4. The duplication part seems fine however, the checks seem to pass and there is no dramatic change in the amount of rows based on merge (none at all, based on my memory). I'd want to check this again, as my old checks may have variable name issues I didn't expect. 

5. Need to be careful with how the files themselves are being saved, so I did not run any "write_metrics" related functions. This also should be pretty carefull assessed because we will need to compile them for output statistics. 

This is also an issue as the actual papermill GPCR bench notebook will need to be updated to account for: 

- Dropping negative and NaN values
- Fixing the CSV reading functions. 

In the future, I should consider some sort of 'don't save' parameter for the metrics reporting. For now, we probably want to just run this in some sort of seperate backup directory and not git add/commit it. Currently, I think it would be best to commit the current state considering how I changed gitignore. 