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

## LIT PCBA 

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
| ADRB2    | 4lde   | 483277     | 110621  | yes   | *KAR - DONE*   |
| ALDH1    | 5l2m   | 228327     | 53846   | yes   | *COS - DONE*   |
| ESR1ago  | 2qzo   | 8341       | 2146    | yes   | *ANT - DONE*   |
| ESR1ant  | 2iog   | 7452       | 1852    | yes   | *ANT - DONE*   |
| FEN1     | 5fv7   | 558263     | 106725  | yes   | *KAR - DONE*   |
| GBA      | 2v3d   | 475989     | 475980  | no    | no    | * no inactives_rdkit.log
| IDH1     | 4umx   | 566613     | 106712  | yes   | *COS - DONE*   |
| KAT2A    | 5mlj   | 540568     | 132219  | yes   | *KAR - DONE*   |
| MAPK1    | 4zzn   | 111544     | 23063   | yes   | *ANT - DONE*   |
| MTORC1   | 4dri   | 41057      | 8003    | yes   | *KAR - DONE*   |
| OPRK1    | 6b73   | 419268     | 92267   | yes   | *COS - DONE*   |
| PKM2     | 3gr4   | 383472     | 110511  | yes   | *ANT - DONE*   |
| PPARG    | 3b1m   | 7751       | 1709    | yes   | *ANT - DONE*   |
| TP53     | 3zme   | 6035       | 1609    | yes   | *ANT - DONE*   |
| VDR      | 3a2j   | 567631     | 107269  | yes   | *KAR - DONE*   |

### Karla: 

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ADRB2/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -JOBNAME ADRB2_inactives_rdkit_ligprep
JobId: karla-0-65f9da6b
```

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/FEN1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:125 -NJOBS 450 -JOBNAME FEN1_inactives_rdkit_ligprep
```

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/KAT2A/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -JOBNAME KAT2A_inactives_rdkit_ligprep
JobId: karla-0-65fa269a
```

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/VDR/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -JOBNAME VDR_inactives_rdkit_ligprep
JobId: karla-0-65fa42bb
```

### Cosmos:

```sh
ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ALDH1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:100 -NJOBS 450 -JOBNAME ALDH1_inactives_rdkit_ligprep
JobId: cosmos-0-65f9e622
```

```sh
ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/IDH1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:100 -NJOBS 450 -JOBNAME IDH1_inactives_rdkit_ligprep
JobId: cosmos-0-65fa0ac8
```

```sh
ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/OPRK1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:100 -NJOBS 450 -JOBNAME OPRK1_inactives_rdkit_ligprep
JobId: cosmos-0-65fa3f73
```

### Tobias:

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/MTORC1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -JOBNAME MTORC1_inactives_rdkit_ligprep
```

### Anton:

```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ESR1_ago/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 80 -JOBNAME ESR1_ago_inactives_rdkit_ligprep
JobId: anton-0-65f9f0e6
```
Note: some weirdness in how the path is handled for ligprep, needs the $SCHRODINGER

```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/ESR1_ant/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 80 -JOBNAME ESR1_ant_inactives_rdkit_ligprep
```
```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/PPARG/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 160  -JOBNAME PPARG_inactives_rdkit_ligprep
JobId: anton-0-65fa0c2a
```

```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/TP53/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 160 -JOBNAME TP53_inactives_rdkit_ligprep
JobId: anton-0-65fa0cb9
```

```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/MAPK1/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 160 -JOBNAME MAPK1_inactives_rdkit_ligprep
JobId: anton-0-65fa27a8
```

```sh
$SCHRODINGER/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/PKM2/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:80 -NJOBS 160 -JOBNAME PKM2_inactives_rdkit_ligprep
JobId: anton-0-65fa3e32
```
#### Quality Control 


- All inactives have re-run, need to rsync the servers. Probably best to sync everything local, make sure the structure is good, and then sync back to server. 

```sh
rsync -avhuzt "karla:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
rsync -avhuzt "cosmos:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
rsync -avhuzt "anton:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
rsync -avhuzt "tobias:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
```
- NOTE: I need to think about updating the rsync script I use because it doesn't do 'u' and 't' options as is. Those options allow for the *most* up to date file to be perserved, which is how I thought rsync worked natively - but it does not. The issue for me is that I might find more footguns I don't expect. However, the basic logic for syncing in this way is probably beneficial. The strong caveat is that you should really be syncing in an intermediate location to confirm that everything is syncing as desired, which thankfully I did here. I also want a --progress option. In general, the CLI version of the script could be really improved. I'd really like to work on this when more critical work is done.

- **NOTE: There is a chance the active ligands failed as well. I should should check those, hopefully since they are small it should be fine.**

In terms of quickly getting results across my files, remember that I have tools like fd, rg, and so on locally. For instance I can do this: 

```sh
fd -e sdf | xargs rg '\$\$\$\$' -c | sort
```

to get something like: 

```sh
ADRB2/actives_rdkit_ligprep.sdf:34
ADRB2/inactives_rdkit.sdf:483268

ALDH1/actives_rdkit_ligprep.sdf:11545
ALDH1/inactives_rdkit.sdf:228324

ESR1ago/actives_rdkit_ligprep.sdf:15
ESR1ago/inactives_rdkit.sdf:8340

ESR1ant/actives_rdkit_ligprep.sdf:135
ESR1ant/inactives_rdkit.sdf:7451

FEN1/actives_rdkit_ligprep.sdf:650
FEN1/inactives_rdkit.sdf:558254

GBA/actives_rdkit.sdf:221
#GBA/actives_rdkit_ligprep.sdf:221
GBA/inactives_rdkit.sdf:475980
#GBA/inactives_rdkit_ligprep.sdf:475980

IDH1/actives_rdkit_ligprep.sdf:61
IDH1/inactives_rdkit.sdf:566604

KAT2A/actives_rdkit_ligprep.sdf:292
KAT2A/inactives_rdkit.sdf:540559

MAPK1/actives_rdkit_ligprep.sdf:403
MAPK1/inactives_rdkit.sdf:111535

MTORC1/actives_rdkit_ligprep.sdf:131
MTORC1/inactives_rdkit.sdf:41056

OPRK1/actives_rdkit_ligprep.sdf:32
OPRK1/inactives_rdkit.sdf:419259

PKM2/actives_rdkit_ligprep.sdf:635
PKM2/inactives_rdkit.sdf:383463

PPARG/actives_rdkit_ligprep.sdf:34
PPARG/inactives_rdkit.sdf:7750

TP53/actives_rdkit_ligprep.sdf:153
TP53/inactives_rdkit.sdf:6034

VDR/actives_rdkit_ligprep.sdf:1442
VDR/inactives_rdkit.sdf:567622
```
- Looks good in terms of inactives 
- Now to check the actives 

```sh
$ fd --glob actives_rdkit.smi
ADRB2/actives_rdkit.smi
ALDH1/actives_rdkit.smi
ESR1ago/actives_rdkit.smi
ESR1ant/actives_rdkit.smi
FEN1/actives_rdkit.smi
GBA/actives_rdkit.smi
IDH1/actives_rdkit.smi
KAT2A/actives_rdkit.smi
MAPK1/actives_rdkit.smi
MTORC1/actives_rdkit.smi
OPRK1/actives_rdkit.smi
PKM2/actives_rdkit.smi
PPARG/actives_rdkit.smi
TP53/actives_rdkit.smi
VDR/actives_rdkit.smi
```

That's 15, good. 

Whare are their `.smi` line counts? 

```sh
fd --glob 'actives_rdkit.smi' | xargs wc -l
```

```sh
35 ADRB2/actives_rdkit.smi
   11546 ALDH1/actives_rdkit.smi
      16 ESR1ago/actives_rdkit.smi
     136 ESR1ant/actives_rdkit.smi
     651 FEN1/actives_rdkit.smi
     222 GBA/actives_rdkit.smi
      62 IDH1/actives_rdkit.smi
     293 KAT2A/actives_rdkit.smi
     404 MAPK1/actives_rdkit.smi
     132 MTORC1/actives_rdkit.smi
      33 OPRK1/actives_rdkit.smi
     636 PKM2/actives_rdkit.smi
      35 PPARG/actives_rdkit.smi
     154 TP53/actives_rdkit.smi
    1443 VDR/actives_rdkit.smi
   15798 total
```

And what do their logs look like?

```sh
fd --glob "actives_rdkit_ligprep.log" | xargs grep "# of processed structures in "actives_rdkit_ligprep.sdf""
```

```sh
ADRB2/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 34
ALDH1/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 11545
ESR1ago/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 15
ESR1ant/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 135
FEN1/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 650
GBA/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 221
IDH1/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 61
KAT2A/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 292
MAPK1/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 403
MTORC1/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 131
OPRK1/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 32
PKM2/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 635
PPARG/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 34
TP53/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 153
VDR/actives_rdkit_ligprep.log:# of processed structures in "actives_rdkit_ligprep.sdf" : 1442
```
It looks like we are all good here. It might be worth to check for "fail"? But its unlikely add much. 

```sh
fd --glob "inactives_rdkit.log" | xargs rg "job\(s\)"
```
```sh
ESR1ago/inactives_rdkit.log:80 of 80 job(s) succeeded; 0 job(s) failed.
TP53/inactives_rdkit.log:159 of 159 job(s) succeeded; 0 job(s) failed.
ESR1ant/inactives_rdkit.log:80 of 80 job(s) succeeded; 0 job(s) failed.
PPARG/inactives_rdkit.log:159 of 159 job(s) succeeded; 0 job(s) failed.
MTORC1/inactives_rdkit.log:447 of 447 job(s) succeeded; 0 job(s) failed.
MAPK1/inactives_rdkit.log:160 of 160 job(s) succeeded; 0 job(s) failed.
ALDH1/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
PKM2/inactives_rdkit.log:160 of 160 job(s) succeeded; 0 job(s) failed.
OPRK1/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
ADRB2/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
IDH1/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
VDR/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
KAT2A/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
FEN1/inactives_rdkit.log:450 of 450 job(s) succeeded; 0 job(s) failed.
```

For the `actives_rdkit_ligprep.log` files, they were always run with `-NJOBS 1` so they don't have an option to fail a single job, just comparing the amount of processed structures, which I have already done. 

It looks like we are good to go. 

Should update the rsf script, really tired of the issues. 

I have fixed the rsf script and I also have a new script to enable persistent connections to a server, just a small QoL improvement. It is in `~/scripts/ssh/establish_connections.sh`. It relies on some changes in my `~/.ssh/config` file as well. 

Now, I think we can run the docking. Here how the files are set up: 

```sh
$ cat FEN1_5fv7_active.sh FEN1_5fv7_inactive.sh
/mnt/data/dk/Schrodinger_adv_2021_1/glide -HOST localhost:1 -NJOBS 1 -OVERWRITE -JOBNAME FEN1_5fv7_active_glide FEN1_5fv7_active.in
/mnt/data/dk/Schrodinger_adv_2021_1/glide -HOST localhost:100 -NJOBS 450 -OVERWRITE -JOBNAME FEN1_5fv7_inactive_glide FEN1_5fv7_inactive.in
```
Both the active and inactive has `.sh` files that will call the `.in` files, let's see how those look: 

```sh
GRIDFILE /mnt/data/dk/work/grids_lit-pcba/FEN1/FEN1_5fv7_structcat-out_complex_prepared__grid.zip
LIGANDFILE /mnt/data/dk/work/lit-pcba/FEN1/actives_rdkit_ligprep.sdf
POSE_OUTTYPE ligandlib_sd
DOCKING_METHOD confgen
PRECISION SP
AMIDE_MODE penal
SAMPLE_RINGS True
EPIK_PENALTIES True
```
```sh
GRIDFILE /mnt/data/dk/work/grids_lit-pcba/FEN1/FEN1_5fv7_structcat-out_complex_prepared__grid.zip
LIGANDFILE /mnt/data/dk/work/lit-pcba/FEN1/inactives_rdkit.sdf
POSE_OUTTYPE ligandlib_sd
DOCKING_METHOD confgen
PRECISION SP
AMIDE_MODE penal
SAMPLE_RINGS True
EPIK_PENALTIES True
```
So based on the GRIDFILE and LIGANDFILE path being the same across the servers (I think it is). We can just run them. 
The tricky thing is to run just the desired `.sh` file I think. I can see via fd that we do have a lot of .sh files, but we only need to run the desired ones. The last thing I am not sure about is if I can them from outside the actual directory. I can try a quick demo with an active file like this: 

```sh
# @karla
/mnt/data/dk/work/grids_lit-pcba
sh FEN1/FEN1_5fv7_active.sh
ERROR
```

No, you have to be in the actual directory because the pathing to the inp file is relative from where you run the `sh` command. Annoying but it is what it is. A future solution is to code the output directory so I can finally be done with this annoying path related issues. For now, I think I just run them manually without coding out logic for it. Big jobs on big servers, rsync at the end. 

| Protein  | PDB_ID | wc -l .smi | rerun glide    |
|----------|--------|------------|----------------|
| ADRB2    | 4lde   | 483277     | *KAR - RUNNING*   |
| ALDH1    | 5l2m   | 228327     | *COS - RUNNING*   |
| ESR1ago  | 2qzo   | 8341       | *ANT - RUNNING*   |
| ESR1ant  | 2iog   | 7452       | *ANT - RUNNING*   |
| FEN1     | 5fv7   | 558263     | *KAR - RUNNING*   |
| GBA      | 2v3d   | 475989     | no      | * no inactives_rdkit.log
| IDH1     | 4umx   | 566613     | *COS - RUNNING*   |
| KAT2A    | 5mlj   | 540568     | *KAR - RUNNING*   |
| MAPK1    | 4zzn   | 111544     | *COS - RUNNING*   |
| MTORC1   | 4dri   | 41057      | *ANT - RUNNING*   |
| OPRK1    | 6b73   | 419268     | *COS - RUNNING*   |
| PKM2     | 3gr4   | 383472     | *ANT - RUNNING*   |
| PPARG    | 3b1m   | 7751       | *ANT - RUNNING*   |
| TP53     | 3zme   | 6035       | *ANT - RUNNING*   |
| VDR      | 3a2j   | 567631     | *KAR - RUNNING*   |

Anton doesn't have the grids_lit_pcba directory, running rsync.

Anton can run the smaller jobs and have Karla/Cosmos run the larger ones. Still, without Tobias, this will take awhile. 

In the meanwhile, we could handle the GPCR-Bench updates and the change we'll need for Pareto ranks on scale (either skipping it temporarily or optimizing it)

Well, EXTREMELY UNFORTUNATELY, anton seems to have died. I either run all those jobs elsewhere, or wait for the server. I updated the lab.

### Anton Re-Runs

| Protein  | PDB_ID | wc -l .smi | rerun glide    |
|----------|--------|------------|----------------|
| ESR1ago  | 2qzo   | 8341       | *ANT - RUNNING kar*   |
| ESR1ant  | 2iog   | 7452       | *ANT - RUNNING kar*   |
| MTORC1   | 4dri   | 41057      | *ANT - RUNNING kar*   |
| PKM2     | 3gr4   | 383472     | *ANT - RUNNING cos*   |
| PPARG    | 3b1m   | 7751       | *ANT - RUNNING kar*   |
| TP53     | 3zme   | 6035       | *ANT - RUNNING kar*   |

Run: PKM2 then MTOR1, then small jobs. 

- Torsion_Strain 

then 

- merge/concat data pipeline for analysis 

## GPCR Bench 

### integrate updated papermill notebook

Need to back up GPCR-bench master first.

```sh
cp -r GPCR-Bench-master/ v2-GPCR-Bench-master
```

file is 

```sh
v2_gpcr_papermill.py
```

We need to convert to a ipynb, inject params, make the papermill runner user it

I have the injected params, need to find the note book that runs it. 

Ran the notebook. Going to inspect some results briefly and then put together a better data analysis pipeline. 

What? All the outputs are in correct. Even though I ran it on the v2 notebook it doesn't have any v2 notebook style data? Oh, I did not actually call v2 where it runs the notebook, I only changed what it inspects. This is now updated. 

It would probably be best to have the inspect/run act on a variable set to the notebook to prevent these issues in the future. 

### Statistics Scripts 

1. Run `v2_GPCR-Bench-master/papermill/csv/combine_data.ipynb` to generate the intermediate and combined data files. NOTE: Be careful to not run it as an unsaved interactive python file, as it will default to running in a directory you might not expect. 

Files generated:  
`csv/concat/*` (strain_enrichment_metrics_concat_<PDB_ID>, strain_log_aucs_concat_<PDB_ID>, strain_roc_metrics_concat_<PDB_ID>)
`csv/combined_data.csv`

`combined_data.csv` has this form: 

| Protein | Strain Energy Cutoff | EF1% | EF5% | deltaEF1% | deltaEF5% | Linear Log10 AUC (x10) | Delta Linear Log10 AUC (x10) | ROC_AUC | Actives | Total Count | deltaAUC |
|---------|---------------------|------|------|-----------|-----------|------------------------|-------------------------------|---------|---------|-------------|----------|

It does not include the 'pareto as scores' style data. 

| Strain Energy Cutoff |
|----------------------|
| No Cutoff |
| 4 |
| 4.5 |
| 5.0 |
| 5.5 |
| 6.0 |
| 7.0 |
| 7.5 |
| 8.0 |
| Top 10 Pareto Ranks |
| Top 20 Pareto Ranks |

The next row would be `No Cutoff` again (you can run the below to check)

```sh
cut -d',' -f2 combined_data.csv | head -n 13 | tail -n +2
# head -n 13 will produce the repeated no cut off line
```

So, we should be able to use this file with the LIT PCBA data as well, as we are unlikely to do the Pareto scores approach. 

Now, we need to recall how the summary statistics are generally found. 

For this, we can run `/v2_GPCR-Bench-master/papermill/combined_data_analysis.ipynb`

This file does *not* do 'pareto as scores' as I think we are moving on from that (I'll check with Francesco). 

This should be a sufficient workflow for LIT PCBA data withone (pretty large) exception. The LIT-PCBA data has an extra underscore tying the `PROTEIN` to the `PDB-ID`. We either need to change how we parse the data itself or we convert that `_` to another character, maybe `-`. This will need to be done on the output stage of the initial csv generation (e.g. the papermill execution), or some sort of `mv` command. Not sure which is best. The extra `mv` command introduces a source of forgetfulness, a refactor is more responsible but will require some sort of in place replace of `title_suffix`. I think that's actually the smartest approach, as it will be operable with the GPCR data (there is no `_` that would convert out). However, the issue is actually earlier, as populating the papermill notebooks will have to be refactored. I think the change will have to be there, because ultimately the `title_suffix` can be replaced via it's very initalization, this will make more sense when I look at the code. I should do that next, as the docking jobs are still running. Maybe I can pull down an rsync version of one server and start to do demo runs. All of them are currently running, so I am not sure which is best. 

I am not even sure I need to pull anything down, because even the "wrong" files will work for this purpose, they just need to exist. Ah, but they need to be converted correctly... I think I did do that, which is where I noticed the mistake in terms of SDF size. Yes, that appears to be true.

See the sdfs:

```sh
$ fd -e sdf
ADRB2/ADRB2_4lde_active_docking_lib_sorted.sdf
ADRB2/ADRB2_4lde_active_glide_lib.sdf
ADRB2/ADRB2_4lde_inactive_docking_lib_sorted.sdf
ADRB2/ADRB2_4lde_inactive_glide_lib.sdf
ALDH1/ALDH1_5l2m_active_glide_lib.sdf
ALDH1/ALDH1_5l2m_inactive_glide_lib.sdf
ESR1ago/ESR1ago_2qzo_active_glide_lib.sdf
ESR1ago/ESR1ago_2qzo_inactive_glide_lib.sdf
FEN1/FEN1_5fv7_active_glide_lib.sdf
FEN1/FEN1_5fv7_inactive_glide_lib.sdf
GBA/GBA_2v3d_active_glide_lib.sdf
GBA/GBA_2v3d_inactive_glide_lib.sdf
IDH1/IDH1_4umx_active_glide_lib.sdf
IDH1/IDH1_4umx_inactive_glide_lib.sdf
KAT2A/KAT2A_5mlj_active_glide_lib.sdf
KAT2A/KAT2A_5mlj_inactive_glide_lib.sdf
MAPK1/MAPK1_4zzn_active_glide_lib.sdf
MAPK1/MAPK1_4zzn_inactive_glide_lib.sdf
MTORC1/MTORC1_4dri_active_glide_lib.sdf
MTORC1/MTORC1_4dri_inactive_glide_lib.sdf
OPRK1/OPRK1_6b73_active_glide_lib.sdf
OPRK1/OPRK1_6b73_inactive_glide_lib.sdf
PKM2/PKM2_3gr4_active_glide_lib.sdf
PKM2/PKM2_3gr4_inactive_glide_lib.sdf
PPARG/PPARG_3b1m_active_glide_lib.sdf
PPARG/PPARG_3b1m_inactive_glide_lib.sdf
TP53/TP53_3zme_active_glide_lib.sdf
TP53/TP53_3zme_inactive_glide_lib.sdf
VDR/VDR_3a2j_active_glide_lib.sdf
VDR/VDR_3a2j_inactive_glide_lib.sdf
```
However, I did not finish running the strain. 

Remember I have a `glide_sort_writer.py` in `/chembl_ligands/grids_lit-pcba`

We should work on making the strain writer now. 

I think rather than trying to think so many steps ahead, it is easier to just make a strain writer that adapts to the current naming conventions. I think refactored at the output of the metric files is most appropriate, as I already have all the names set from the docking runs. I could also re-write the glide sort utility to use the hyphen, and then it is perserved. Let's try it like that, since we will need to re-run all the glide sort utility as well. 

I should also get a workflow to get rid of "v2" I think I can just `cp` or `rsync` over my results, but we can do that later. 

`glide_sort_writer.py` has been updated to write PROTEIN-PDBID files, which hopefully reduce complexity moving forward. However, I did replace `.sdfgz` with `_sorted.sdf` which may have been a bad idea. 

It is best to now write the strain script. A little worried about the anti-overwriting functionality right now in my torsion scripts. Probably want to either find a careful `rm` command for prototyping. 

I have `strain_writer.py` at a point where it is about ready to run the commands, I think I just need to integrate something like this: 

```py
import subprocess

# Assuming `commands` is your list of commands to run in parallel
# For example: commands = ['conda run -n myenv python path/to/script.py', 'conda run -n myenv python path/to/script.py']
commands = [...]
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
stdout, stderr = process.communicate(input=commands_string)

# Optionally, handle the output and errors
print(stdout)
if stderr:
    print(f"Errors:\n{stderr}")

# Example for saving stdout and stderr to files
with open('parallel_strain_writer.log', 'w') as f_out, open('parallel_errors.txt', 'w') as f_err:
    f_out.write(stdout)
    if stderr:
        f_err.write(stderr)
```
Best to run this from home/server to avoid interruptions. 

Karla has `glide_sort` on path, should be good to run directly. Need to rsync local from server and then onto karla.

Ran karla -> local. Once cosmos is done I can do that. I want to double check GBA though. 

From my notes, I put everything related to ligand preparation: 

```sh
rsync -avhuzt "karla:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
rsync -avhuzt "cosmos:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
rsync -avhuzt "anton:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
rsync -avhuzt "tobias:'/mnt/data/dk/work/lit-pcba/'" "/Users/lkv206/work/lit-pcba/"
```

If go to that directory and check: 

```sh
/Users/lkv206/work/lit-pcba/GBA
$ grep processed inactives_rdkit_ligprep.log | tail -n 1
# of processed structures in "inactives_rdkit_ligprep.sdf" : 475980
```
Note that while this is `GBA`, it has the 'wrong' file ending. So the docking file should be checked as well. 

```sh
/mnt/data/dk/Schrodinger_adv_2021_1/glide -HOST localhost:100 -NJOBS 450 -OVERWRITE -JOBNAME GBA_2v3d_inactive_glide GBA_2v3d_inactive.in
```

```sh
$ cat GBA_2v3d_inactive.in
GRIDFILE /mnt/data/dk/work/grids_lit-pcba/GBA/GBA_2v3d_structcat-out_complex_prepared__grid.zip
LIGANDFILE /mnt/data/dk/work/lit-pcba/GBA/inactives_rdkit.sdf
POSE_OUTTYPE ligandlib_sd
DOCKING_METHOD confgen
PRECISION SP
AMIDE_MODE penal
SAMPLE_RINGS True
EPIK_PENALTIES True
```
Huh, that is saying it ran on inactive_rdkit.sdf, also I can it ran on `tobias`. According to rsync I did pull down tobias, so let me check for `/mnt/data/dk/work/lit-pcba/GBA/inactives_rdkit.sdf`.

Ah, right, no log file for the ligprep... worth it to rerun? Maybe I can find it on Tobias. 

No I can't find it. Honestly, I think the re-run and re-dock is better for my sanity? But it looks SO correct:

```sh
@karla 
/mnt/data/dk/work/lit-pcba/GBA
❯ wc -l inactives_rdkit.smi
475989 inactives_rdkit.smi
❯ grep "\$\$\$\$" inactives_rdkit.sdf | wc -l
475980
```

I think it is fine to run as is. If there is something wrong with the GBA data, I'll look at it again. I can compare stuff between the .smi and .sdf to make sure. Do I have the commands to run LigPrep and Glide easily? Glide it seems like for sure. But b/c there is no log for the LigPrep, I will have to find and re-run the command. 

There is just no way I "accidentally" made the exactly correct file after starting with an incorrect file. Let's just sync for now and we can update later if we need to do. 

Need cosmos sync to local. 

Ah, cosmos is still running... 

I guess we may as well do the LigPrep then, it will likely be on par with whatever cosmos is doing. Now I have these records to tell me that I AM RERUNNING LIGPREP FOR GBA ON KARLA. If the results are the same, chance is even less. 

```sh
/mnt/data/dk/Schrodinger/ligprep -inp /mnt/data/dk/scripts/job_writer/ligprep.inp -ismi /mnt/data/dk/work/lit-pcba/GBA/inactives_rdkit.smi -osd inactives_rdkit.sdf -HOST localhost:150 -NJOBS 450 -JOBNAME GBA_inactives_rdkit_ligprep
```

Okay, I launched it. 

Strangely, karla's `/mnt/data/dk/work/grids_lit-pcba` doesn't have `.sh` files for ONLY `GBA`. 

Just gonna confirm I have that locally and then run in reverse to get it there:

```sh
rsync -avhuzt --progress "karla:'/mnt/data/dk/work/grids_lit-pcba/'" "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/"
Proceed with the transfer? (y/n): y
receiving file list ...
703 files to consider

sent 16 bytes  received 19.02K bytes  12.69K bytes/sec
total size is 10.96G  speedup is 575961.98
```

All good. So: 

```sh
rsync -avhuzt --progress "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/" "karla:'/mnt/data/dk/work/grids_lit-pcba/'"
```

```sh
building file list ...
900 files to consider
adrb2_duplicate_check.py
67.86K 100%   33.47MB/s    0:00:00 (xfer#1, to-check=898/900)
```
I should build a logging mechanism. Regardless, this is good. 

So we'll get files up there. 

I think my rsync command isn't ideal here because all of the ownership is listed as belonging to lkv206. I would think it wouldn't actually be a very large file transfer because karla/tobias/etc would have more updated copies of these files? No, unfortunately rsync isn't very good about perserving ownership and because there is a time difference between me and copenhagen, it is going to do the reupload again. Then on karla, it will be owned by 'dk' and lkv206, so a transfer locally will again be painful unless I specify some sort of time argument. In theory, when it is 9AM est, it is 3PM in copenhagen where the servers are, so copenhagen is "newer" by default. Why would it trigger a re-transfer then in my case? It seems it doesn't work like how I expect. Maybe there is some sort of checksum thing I can do, or maybe like git can do this somehow? I have no idea. If it becomes a big issue we can handle this. It's either this or git LFS, and I hate git lfs. 

It seems checksum is best 

```sh
rsync -avhuzt --checksum --progress "karla:'/mnt/data/dk/work/grids_lit-pcba/'" "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/"
```

I will integrate this later. 

GBA glide jobs started on Karla.

Cosmos rsync to local started. 

Once both are done, (1) all to local, (2) push to karla, (3) glide sort, (4) strain, (5) fix notebook, execute notebook, (6) summary stats. 

1. All to local:

 Karla -> local 

```sh
rsync -avhuzt --checksum --progress "karla:'/mnt/data/dk/work/grids_lit-pcba/'" "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/"
```

Cosmos -> local 

I need to wait until the Karla -> Local is done to run the below, hard to interpert the --dry-run output. 

```sh
rsync -avhuzt --checksum --progress "cosmos:'/mnt/data/dk/work/grids_lit-pcba/'" "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/"
```
Both done. 

2. Push to Karla

```sh
rsync -avhuzt --checksum --progress "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/" "karla:'/mnt/data/dk/work/grids_lit-pcba/'"
```
Really wish the transfer speed was faster. 

3. glide sort:

I probably run this without much issue using karla's `data_viz` conda env

Oh, I think there is an issue, we are not getting the updated `protein-pdbid`, this is because the amount of underscores `_` is different. 

Need to refactor this approach to just the base name somehow. 

Very satisfying, figured that out almost completely by myself. Used the os library to handle pathing instead of naively counting underscores and seperating them. 

Use git to commit the change with a proper message, rsync to karla. 

Now to delete unwanted sdf files, to prevent clutter. 

```sh
fd -e sdf | xargs rm
```
I should run this locally as well to clean the directory. 

It would be nice to have a sort of `backup.sh` script that automatically backs up my directory to some location, but that's a later project. 

You should use `--exec` when possible, e.g:

```sh
fd -e sdf --exec rm {}
```
Strangely, the output files didn't filter any poses out? Weird but not impossible, will simply be careful when examining papermill output.

4. strain_writer:

I should also update `strain_writer.py` to use more portable logic? But it actually seems okay right now.

The bigger issue is that I want to have analytics_env on karla as I will otherwise have portability issues. 

Fixed all the env/execution related issues found so far, executed the script. It is not very fast, waiting on completion, then need a karla -> local rsync. Then ssh papermill set up. 

Local sync: 

```sh
rsync -avhuzt --checksum --progress "karla:'/mnt/data/dk/work/grids_lit-pcba/'" "/Users/lkv206/work/to_do_projects/chembl_ligands/grids_lit-pcba/"
```

5. papermill 

The way papermill works (at least my set up) is having a "runner", "template" and "output". The runner acts on template and injects the parameters (e.g. the variables that are going to be changing) into the template. The runner then submits the template for execution, the template has logic for defining the output ipynb. The enrichment analysis code is run in the output files. 

Here is how it works in GPCR-bench-master/

Runner: `GPCR-bench-master/papermill.ipynb`

Template: `GPCR-bench-master/gpcr_papermill.ipynb`

Output: `GPCR-bench-master/papermill/notebooks/gpcr_papermill_output_{title_suffix}.ipynb`

We can use a similar structure for LIT PCBA, but changing the obvious project variable names, and then defining appropriate logic for the PROTEIN-PDBID notation we now have.

Note: Had to do some debugging for the file paths, but have an interactive python window of `litpcba_papermill.py` running. 

New set up for `grids_lit-pcba`

Runner: `grids_lit-pcba/papermill.ipynb`

- Need to update logic. 

Template: `grids_lit-pcba/litpcba_papermill.ipynb`

Output: what you expect 

Looooooots of debugging. Code seems to be running locally, will rsync to Karla. If it finishes a notebook, I will run all of it on Karla too.

Running on karla, not sure if it's gonna quit on nohup? 

Execution on Karla is very slow... Not sure what's going on there. 

6. Summary Statistics 

I am really suspicious about just the paperill output because performance is really low. Did I prepare the grids correctly via the automated workflow? How does my performance compare to the paper? For now, let's get the summary stats going. 

Step 1 (from GPCR-bench) was `papermill/csv/combine_data.ipynb` to make `papermill/csv/combined_data.csv`

Step 2: `GPCR-Bench-master/papermill/combined_data_analysis.py`

Ran this all for `grids_lit-pcba` as well (same naming convention). Had to delete erroneous data that was generated prior to updated the title_suffix convention to use a hyphen (`-`). Should delete on Karla as well, just in case. 

@Karla:

```sh
fd --glob "*ADRB2*
```

```sh
csv/strain_enrichment_metrics_ADRB2-4lde.csv
csv/strain_enrichment_metrics_pareto_ADRB2-4lde.csv
csv/strain_log_aucs_ADRB2-4lde.csv
csv/strain_log_aucs_pareto_ADRB2-4lde.csv
csv/strain_roc_metrics_ADRB2-4lde.csv
csv/strain_roc_metrics_pareto_ADRB2-4lde.csv
notebooks/ADRB2-4lde_litpcba_papermill.ipynb
```

Looks good, deletes are done. 

I should write up the data or simply make the plots look nicer when possible. 

However, data seems to imply a cut off of 6 or so is good. This is pretty similar the GPCR-bench data as well. 

Next steps: 

Double check docking
Cheminformatics type analysis
Deep Docking 

All are valid workflows. Only docking double check will be done here. 

## Docking Check:

Unfortunately, need to unzip the files to load the receptor. 

I loaded all the required files into a Maestro project file:

`grids_lit-pcba/docking_complexes/all_complexes_actives.prj`

It is fighting with me over permissions to export it as a sdfgz. Need to be connected to the VPN. 

The complexes are exported now. 

`grids_lit-pcba/docking_complexes/all_complexes_actives.sdfgz`

I guess the most practical analysis is visual comparison. I don't think it's going to be very straight forward. 

Co-Crystal Binding Pocket: Do they mostly overlap with top score active? 

ADRB2-4lde: Yes

ALDH1-5l2m: Yes
- I don't know where the mix up happened, but *5l2m* is the correct PDB for this complex. THe structure I used is correct, the binding site overlaps. Just becareful later.

ESR1ago-2qzo: Yes
- Note this structure is a dimer with some disconnected loops that look like additional ligands, but I still have the correct assignment. 

ESR1ant-2iog: Yes
- Note that looks like a hard pocket, active would need to form a "V" shape (think o/m subsitution, cis-alkene) to fit.

FEN1-5fv7: Yes
- While the grid is definitely in the correct location, the overlap of the top scoring active and cocrystal is not great. The scores of the actives in general are pretty bad. However, it seems very unlikely that the grid itself was prepared incorrectly. 

GBA-2v3d: Yes
- Also a dimer. Clearly a challenging pocket I would say. We are however docking to the correct location, very strong overlap. 

IDH1-4umx: Yes 

KAT2A-5mlj: Yes 
- Dimer 

MAPK1-4zzn: Yes
- This exercise is really leading me to believe in the "unhappy water hypothesis", I feel like a lot of these systems are not accounting for them, but I can't be sure they are relevant. 

MTORC1-4dri: Yes
- Yeah it's definitely the right spot but the binding pocket FK506! It's the rapamyacin thing that Schreiber is famous for. The brief actives I looked at are not marcocycles like the cocrystal ligand is. I have to imagine the strain here is huge, and it would be in the actives as well. 

OPRK1-6b73: Yes 

I am also unsurprised here. Opiates are weird. Though I guess fentanyl exists so you can capture the pharmacophore without really needing the bicyclic rings. 

PKM2-3gr4: Yes
- Yeah these waters are perflectly happy being displaced. Weird structure though. The ligands remind me of my project. 

PPARG-3b1m: Yes

TP53-3zm3: Yes

VDR-3a2j

All of them passed. Yes!

Now I want to clean up the directory to make this project easier to read. 

Some deletions are trivial, others maybe less so. I am going to cp v2_GPCR bench to GPCR bench.

I actually need to mv. Do I need the old GPCR-bench-master? probably not. The only real advantage is that I have the ParetoRanks data (at best). I have my back up just in case anyway so let's just move it. 




## Summary and Notes:


