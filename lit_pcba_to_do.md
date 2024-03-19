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

# Strain Metrics Cont. 

I replicated the GPCR_bench-master/ directory structure for grids_lit-pcba, so that I could save my files in a similar way, rather than refactoring the entire saving logic to instead using some sort of "directory" variable as that set up is very tedious. 

I am now hopefully able to progress through the rest of the interactive python file, once that is complete, I can return to parallel execution of the sort/strain calculations with my prior logic. I am a little worried about memory allocation when running the strain calculations.

# Strain Metrics Intermission

Testing the final section, "Pareto Ranks as Scores". This calculation will take awhile, so we can focus on other tasks we need to complete right now: 

# # LIT PCBA 

- `glide_sort`

See: `grids_lit-pcba/glide_sort_writer.py`

- - NOTE: Naming convention changed for the output files. I am simply replacing `.sdfgz` for `.sdf` between input and output, I forgot to add the `_sorted` prefix I originally planned for. This is fine just be aware of it.
- - NOTE: While some files (e.g. ADRB2) have quite a large amount of inactives (>100K), other targets are *considerably* smaller. See this example: 

```
Glide Sort
------------------------------------------------------------------------------
REPORT OF BEST 7759 POSES
The sorted ligand structures were written to the file:
/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/MTORC1/MTORC1_4dri_inactive_glide_lib.sdf

```

This is much smaller than what is mentioned on their website for MTORC1: 

| AID   | Set    | Target | Ligands | Actives | Inactives | PDB templates |
|-------|--------|--------|---------|---------|-----------|---------------|
| 493208| MTORC1 | Mechanistic target of rapamycin | Inhibitors | 97 | 32972 | 11 |

It is very much possible to lose a fair amount of ligands to the binding pocket itself (I have seen this in my own research, such as a Deep Docking run of GPR3, where I was docking less than 20% of my intended molecules). But this is a pretty sharp decrease. I think its worth investigating just to double check. 

However, when investigating the actual job log. I see something pretty strange: 

`Some jobs failed: 7765 of 8003 ligands done.`

The failure isn't that weird (I don't think), but why are there only 8K ligands? 

The command is: 

`Commandline    : /mnt/data/dk/Schrodinger_adv_2021_1/glide -HOST localhost:100 -NJOBS 450 -OVERWRITE -JOBNAME MTORC1_4dri_inactive_glide MTORC1_4dri_inactive.in`

If we check the in file: 

```sh
$ cat MTORC1_4dri_inactive.in
GRIDFILE /mnt/data/dk/work/grids_lit-pcba/MTORC1/MTORC1_4dri_structcat-out_complex_prepared__grid.zip
LIGANDFILE /mnt/data/dk/work/lit-pcba/MTORC1/inactives_rdkit.sdf
POSE_OUTTYPE ligandlib_sd
DOCKING_METHOD confgen
PRECISION SP
AMIDE_MODE penal
SAMPLE_RINGS True
EPIK_PENALTIES True
```
Now let's go check `inactives_rdkit.sdf` on Tobias

Examining the log (`inactives_rdkit.log`) for the LigPrep job returns: 

`Commandline    : /mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/MTORC1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -WAIT -JOBNAME MTORC1_inactives_rdkit_ligprep`

Taking a line count of our input: `/mnt/data/dk/work/lit-pcba/MTORC1/inactives_rdkit.smi`

```sh
$ wc -l inactives_rdkit.smi
    41057 inactives_rdkit.smi
```
So 41K, which is reasonable given the stereochemical enumeration. The failed ligands seem to be because of an overuse of licenses during LigPrep. We can see that in `inactives_rdkit.log`, towards the end of the file: 

```
--------------------------------------------------------------------------------
JobId          : karla-0-65b870ce
Name           : inactives_rdkit-444
Program        : LigPrep
Version        : 2022-2 build 128
Host           : karla.sund.root.ku.dk
Dir            : /mnt/data/dk/scratch/dk/inactives_rdkit.9
HostEntry      : localhost
JobHost        : karla.sund.root.ku.dk
JobDir         : /mnt/data/dk/scratch/dk/inactives_rdkit-444.4
JobSchrodinger : /mnt/data/dk/Schrodinger
Commandline    : /mnt/data/dk/Schrodinger/ligprep -ma 500 -bff 16 -i 0 -s 1 -nd -nt -orig_file /mnt/data/dk/work/lit-pcba/MTORC1/inactives_rdkit.smi -orig_file_index 40757 -ismi in_inactives_rdkit-444.smi -osd inactives_rdkit-444.sdf -HOST localhost
StartTime      : 2024-01-30-04:45:19
--------------------------------------------------------------------------------
guard: WARNING Unable to checkout 1 LIGPREP_MAIN v62 license(s). Error is -4: Licensed number of users already reached. Tried 3 time(s) in 30 seconds.
guard: FATAL -4: Licensed number of users already reached.
ERROR: ligprep3_initialize() failed
--------------------------------------------------------------------------------
SUBJOB: inactives_rdkit-445
```
I am guessing this was a general issue. Need to update the team and fix this. 

I updated the team and I will need to start re-running both the LigPrep and the follow up docking job. I wonder if I will experience issues with respect to an --OVERWRITE flag, I sincerely hope not. 

For now, I need a table of some sort where I can catalogue the execution of these ligprep jobs. I have moved files around to allow them to run on Tobias in addition to Karla by getting the .inp file in the identical location. I need to be able to rsync the servers to one another. Probably use cosmos or anton as well if I can. First, I need the rsync tool. 

NOTE: We need to drop the `-WAIT` flag for now I think, otherwise we need to run some sort of `nohup &` process.
NOTE: `--OVERWRITE` isn't a valid command anyway, so I think we are safe on that end, at least. 

Some jobs ran just fine it seems, though it is worthwhile to double check. I am getting a local rsync copy of `lit-pcba/` from tobias. I'll analyze it locally. My rsync server script is not good, but the local one is and I can work locally much faster, so I rather do it all here and push to the servers.

I am going to rsync from Karla as well, it should be the most updated files. 

Surprising result from Karla, might just be because I am rerunning MTORC, but either way: 

```
rsync -avhz "karla:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
Proceed with the transfer? (y/n): y
receiving file list ... done
GBA/
MTORC1/
MTORC1/inactives_rdkit-dropped-indices.txt
MTORC1/inactives_rdkit-dropped.smi
MTORC1/inactives_rdkit.log
MTORC1/inactives_rdkit.sdf

sent 11.89K bytes  received 6.42M bytes  858.24K bytes/sec
total size is 7.34G  speedup is 1140.84
```
Time to start the review...

| Protein  | PDB_ID | wc -l .smi | sdf.log | fails | rerun |
|----------|--------|------------|---------|-------|-------|
| ADRB2    | 4lde   | 483277     | 110621  | yes   | tbd   |
| ALDH1    | 5l2m   | 228327     | 53846   | yes   | tbd   |
| ESR1ago  | 2qzo   | 8341       | 2146    | yes   | tbd   |
| ESR1ant  | 2iog   | 7452       | 1852    | yes   | tbd   |
| FEN1     | 5fv7   | 558263     | 106725  | yes   | tbd   |
| GBA      | 2v3d   | 475989     | 475980  | no    | no    | * no inactives_rdkit.log
| IDH1     | 4umx   | 566613     | 106712  | yes   | tbd   |
| KAT2A    | 5mlj   | 540568     | 132219  | yes   | tbd   |
| MAPK1    | 4zzn   | 111544     | 23063   | yes   | tbd   |
| MTORC1   | 4dri   | 41057      | 8003    | yes   | tbd   |
| OPRK1    | 6b73   | 419268     | 92267   | yes   | tbd   |
| PKM2     | 3gr4   | 383472     | 110511  | yes   | tbd   |
| PPARG    | 3b1m   | 7751       | 1709    | yes   | tbd   |
| TP53     | 3zme   | 6035       | 1609    | yes   | tbd   |
| VDR      | 3a2j   | 567631     | 107269  | yes   | tbd   |

*GBA is likely fine, I think the lack of log was because it was a system I was testing with in another dir and then cp or mv the relevant SDF files over.

Well, unfortunately, this means re-running most of these through both LigPrep and Glide Docking. To save time, need to put split.seperate them by server to get as many parallel CPU jobs as possible. 


| Protein  | PDB_ID | wc -l .smi | sdf.log | fails | rerun |
|----------|--------|------------|---------|-------|-------|
| ADRB2    | 4lde   | 483277     | 110621  | yes   | *KAR*   |
| ALDH1    | 5l2m   | 228327     | 53846   | yes   | *COS*   |
| ESR1ago  | 2qzo   | 8341       | 2146    | yes   | *ANT*   |
| ESR1ant  | 2iog   | 7452       | 1852    | yes   | tbd   |
| FEN1     | 5fv7   | 558263     | 106725  | yes   | tbd   |
| GBA      | 2v3d   | 475989     | 475980  | no    | no    | * no inactives_rdkit.log
| IDH1     | 4umx   | 566613     | 106712  | yes   | tbd   |
| KAT2A    | 5mlj   | 540568     | 132219  | yes   | tbd   |
| MAPK1    | 4zzn   | 111544     | 23063   | yes   | tbd   |
| MTORC1   | 4dri   | 41057      | 8003    | yes   | *TOB*   |
| OPRK1    | 6b73   | 419268     | 92267   | yes   | tbd   |
| PKM2     | 3gr4   | 383472     | 110511  | yes   | tbd   |
| PPARG    | 3b1m   | 7751       | 1709    | yes   | tbd   |
| TP53     | 3zme   | 6035       | 1609    | yes   | tbd   |
| VDR      | 3a2j   | 567631     | 107269  | yes   | tbd   |

### Karla: 

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ADRB2/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -JOBNAME ADRB2_inactives_rdkit_ligprep
JobId: karla-0-65f9da6b
```

### Cosmos:

```sh
ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ALDH1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:100 -NJOBS 450 -JOBNAME ALDH1_inactives_rdkit_ligprep
JobId: cosmos-0-65f9e622
```

### Tobias:

MTORC job

### Anton:

```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ESR1_ago/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 80 -JOBNAME ESR1_ago_inactives_rdkit_ligprep
JobId: anton-0-65f9f0e6
```
Note: some weirdness in how the path is handled for ligprep, needs the $SCHRODINGER

- NOTE: There is a chance the active ligands failed as well. I should should check those, hopefully since they are small it should be fine.


- Torsion_Strain 

then 

- merge/concat data pipeline for analysis 

## GPCR Bench 

- integrate updated papermill notebook 
- rerun data analysis 

