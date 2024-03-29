# Deep Docking Set Up 

## Deep Docking Architecture

Before beginning any change to Deep Docking, it is necessary to have a better understanding of how it works in terms of both general architecture and code. I don't know the best way to approach this without excess redundancy, but I think following the code and `phase_X` workflow is the best way to make a cohesive summary. 

### Phase 1

`phase_1.sh` begins after preparation of the underlying chemical library you are screening against and the set up of your project directory/logs file; it is skipped here for brevity. The log file is presented, as it will be referenced throughout: 

`logs.txt`:

```sh 
/home/DeepDocking/projects                               #Path to the project folder
protein_test_automated                                   #Name of project folder
/home/DeepDocking/docking_grid/glide_grid.zip            #Location of the docking grid
/home/DeepDocking/library_prepared_fp                    #Location of the fingerprint library
/home/DeepDocking/library_prepared                       #Location of the SMILES library
Glide                                                    #Name of the docking program (either Glide or FRED)
24                                                       #Number of models to train (different combinations of hyperparameters)
1000000                                                  #Size of validation and test sets
~/DeepDocking/DD_protocolscripts_1/glide_template.in     #Location of the Glide docking template script (leave blank if FRED is used)
```

 Let's begin by describing `phase_1.sh`. 

 It is run like so, with my annotations in `()` for clarity: 

 ```sh
 sbatch --cpus-per-task n_cpus_per_node phase_1.sh ($1) current_iteration ($2) n_cpus_per_node  ($3) path_project ($4) project ($5) training_sample_size ($6) conda_env
 ``` 

 The script is below: 

 ```sh
 #!/bin/bash
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60       #change this to match the number of processors you want to use
#SBATCH --mem=0
#SBATCH --job-name=phase_1

source ~/.bashrc
conda activate $6

start=`date +%s`

file_path=`sed -n '1p' $3/$4/logs.txt`
protein=`sed -n '2p' $3/$4/logs.txt`
n_mol=`sed -n '8p' $3/$4/logs.txt`

pr_it=$(($1-1)) 

t_cpu=$2

mol_to_dock=$5

if [ $1 == 1 ]
then 
	to_d=$((n_mol+n_mol+mol_to_dock))
else
	to_d=$mol_to_dock
fi
 
echo $to_d
echo $t_cpu

python jobid_writer.py -pt $protein -fp $file_path -n_it $1 -jid $SLURM_JOB_NAME -jn $SLURM_JOB_NAME.txt

morgan_directory=`sed -n '4p' $3/$4/logs.txt`
smile_directory=`sed -n '5p' $3/$4/logs.txt`
sdf_directory=`sed -n '6p' $3/$4/logs.txt`

if [ $1 == 1 ];then pred_directory=$morgan_directory;else pred_directory=$file_path/$protein/iteration_$pr_it/morgan_1024_predictions;fi

python scripts_1/molecular_file_count_updated.py -pt $protein -it $1 -cdd $pred_directory -t_pos $t_cpu -t_samp $to_d
python scripts_1/sampling.py -pt $protein -fp $file_path -it $1 -dd $pred_directory -t_pos $t_cpu -tr_sz $mol_to_dock -vl_sz $n_mol
python scripts_1/sanity_check.py -pt $protein -fp $file_path -it $1
python scripts_1/extracting_morgan.py -pt $protein -fp $file_path -it $1 -md $morgan_directory -t_pos $t_cpu
python scripts_1/extracting_smiles.py -pt $protein -fp $file_path -it $1 -smd $smile_directory -t_pos $t_cpu

end=`date +%s`
runtime=$((end-start))
echo $runtime
```

Line by line, the scripts starts by using information about the project to call arguments for variables we will need from the log file (`($3) path_project ($4) project)`). 

Most of what is happening is straight forward, though to be more clear `pr_it` is referring to the previous iteration and is used when `$1` (current iteration) is not equal to 1. 

The first piece of logic sets the `to_d` variable, presumably "to dock". If we are in iteration 1, then we take `$mol_to_dock` (set when we pass `phase_1.sh` with the required arguments) as well as `n_mol`, which is taken from the `logs.txt` file. If we are *not* in iteration 1, the amount of molecules we dock (`to_d`) is just the `$mol_to_dock`, meaning that the `logs.txt` file is used to set the amount of compounds we will be docking in iterations > 1 (you can of course edit the logs.txt file if you wish).

We echo some information to console and then call: 

```sh
python jobid_writer.py -pt $protein -fp $file_path -n_it $1 -jid $SLURM_JOB_NAME -jn $SLURM_JOB_NAME.txt
```

From what I understand, despite not explicitly setting `$SLURM_JOB_NAME`, we can use that variable due to the SLURM header at the start. 

I will not cover `jobid_writer.py` in much detail here except to say I believe it just handles creating the job directories we will write to. It is pretty straight forward. 

Then, we read the logs.txt file in order to get our `morgan_directory`, `smile_directory`, and `sdf_directory`. Though it does appear that `sdf_directory` is really just setting Glide or FRED. 

Then we hae logic regarding what our prediction directory is (`pred_directory`) is. Iteration 1 will use the `morgan_directory` that is found from the `logs.txt` file. Other wise, it will use the previous iteration (`pr_it`) `morgan_1024_predictions/` directory.

Now, we can run the five python scripts in sequential order: 

```sh
python scripts_1/molecular_file_count_updated.py -pt $protein -it $1 -cdd $pred_directory -t_pos $t_cpu -t_samp $to_d
python scripts_1/sampling.py -pt $protein -fp $file_path -it $1 -dd $pred_directory -t_pos $t_cpu -tr_sz $mol_to_dock -vl_sz $n_mol
python scripts_1/sanity_check.py -pt $protein -fp $file_path -it $1
python scripts_1/extracting_morgan.py -pt $protein -fp $file_path -it $1 -md $morgan_directory -t_pos $t_cpu
python scripts_1/extracting_smiles.py -pt $protein -fp $file_path -it $1 -smd $smile_directory -t_pos $t_cpu
```

Let's start with:

```sh
python scripts_1/molecular_file_count_updated.py -pt $protein -it $1 -cdd $pred_directory -t_pos $t_cpu -t_samp $to_d
```


```py
from multiprocessing import Pool
from contextlib import closing
import pandas as pd
import numpy as np
import argparse
import glob
import time

try:
    import __builtin__
except ImportError:
    # Python 3
    import builtins as __builtin__

# For debugging purposes only:
def print(*args, **kwargs):
    __builtin__.print('\t molecular_file_count_updated: ', end="")
    return __builtin__.print(*args, **kwargs)

parser = argparse.ArgumentParser()
parser.add_argument('-pt','--project_name',required=True,help='Name of the DD project')
parser.add_argument('-it','--n_iteration',required=True,help='Number of current DD iteration')
parser.add_argument('-cdd','--data_directory',required=True,help='Path to directory contaning the remaining molecules of the database ')
parser.add_argument('-t_pos','--tot_process',required=True,help='Number of CPUs to use for multiprocessing')
parser.add_argument('-t_samp','--tot_sampling',required=True,help='Total number of molecules to sample in the current iteration; for first iteration, consider training, validation and test sets, for others only training')
io_args = parser.parse_args()


protein = io_args.project_name
n_it = int(io_args.n_iteration)
data_directory = io_args.data_directory
tot_process = int(io_args.tot_process)
tot_sampling = int(io_args.tot_sampling)

print("Parsed Args:")
print(" - Iteration:", n_it)
print(" - Data Directory:", data_directory)
print(" - Sampling Size:", tot_sampling)


def write_mol_count_list(file_name,mol_count_list):
    with open(file_name,'w') as ref:
        for ct,file_name in mol_count_list:
            ref.write(str(ct)+","+file_name.split('/')[-1])
            ref.write("\n")


def molecule_count(file_name):
    temp = 0
    with open(file_name,'r') as ref:
        ref.readline()
        for line in ref:
            temp+=1
    return temp, file_name


if __name__=='__main__':
    files = []
    for f in glob.glob(data_directory+'/*.txt'):
        files.append(f)
    print("Number Of Files:", len(files))

    t=time.time()
    print("Reading Files...")
    with closing(Pool(np.min([tot_process,len(files)]))) as pool:
        rt = pool.map(molecule_count,files)
    print("Done Reading Finals - Time Taken", time.time()-t)

    print("Saving File Count...")
    write_mol_count_list(data_directory+'/Mol_ct_file_%s.csv'%protein,rt)
    mol_ct = pd.read_csv(data_directory+'/Mol_ct_file_%s.csv'%protein,header=None)
    mol_ct.columns = ['Number_of_Molecules','file_name']
    Total_sampling = tot_sampling
    Total_mols_available = np.sum(mol_ct.Number_of_Molecules)
    mol_ct['Sample_for_million'] = [int(Total_sampling/Total_mols_available*elem) for elem in mol_ct.Number_of_Molecules]
    mol_ct.to_csv(data_directory+'/Mol_ct_file_updated_%s.csv'%protein,sep=',',index=False)
    print("Done - Time Taken", time.time()-t)
```

#### Debugging Print Function
It overrides the built-in `print` function for debugging, prefixing output with "\t molecular_file_count_updated: " to differentiate script output from other messages.

#### Argument Parsing
The five required command-line arguments:
- **project_name** (`-pt`): Name of the DD project.
- **n_iteration** (`-it`): The current iteration number of the DD process.
- **data_directory** (`-cdd`): The path to the directory containing molecule database files.
- **tot_process** (`-t_pos`): The number of CPUs to use for multiprocessing.
- **tot_sampling** (`-t_samp`): The total number of molecules to sample in the current iteration.

All are set from the information we pass in `phase 1.sh`. Importantly, `t_samp` is found via `mols_to_dock`.

#### Main Functionality
1. **File Discovery**: Uses `glob.glob` to find all `.txt` files in the specified data directory as the `files` list.  

2. **Parallel Molecule Counting**:

   - Initializes a multiprocessing pool with the lesser of the number of CPUs requested or the number of `files` to process. This is handled by the `np.min` function: `np.min([tot_process,len(files)])`.

   - The use of `closing`  (from `import contextlib`) ensures that the pool is closed after the `with` block is executed. This is important because the pool is not garbage collected by default (I believe an error could stop execution).

   - Each file is processed in parallel by the `molecule_count` function, which counts the non-header lines in the file. After opening the file in read-only mode, the use of `ref.readline()` reads the first line, but does nothing with it, effectively skipping it (maybe easier to say "consuming it"). Then, we can iterate through the lines (python's default behavior for file objects) and count them. The `temp` variable is used to store the count, and is returned along with the file name. 
   
   - While the syntax is unfamiliar to me, the use of pool.map(function, iterable) is used to map the execution of a function (`molecule_count`) to an iterable (`files`). From reading into it, this has to do with functions in python being first-class objects, and the `multiprocessing` module being able to use this to distribute work. You could replace `molecule_count` with something like `print` and see that it is called for each file in `files`. This is stored as a tuple in the `rt` list. For example `rt[0]` is the tuple `(count, file_name)` for the first file in `files`, `rt[0][0]` is the count (read it as "at list 0, element 0 within the tuple at list 0"), and `rt[0][1]` is the file name (read it as "at list 0, element 1 within the tuple at list 0"). 

   - Thus, we get a list of tuples, where each tuple is the line count and file name for a given file. 

3. **Result Compilation and Saving**:
   - The molecule counts and file names are written to a CSV file (`data_directory + '/Mol_ct_file_' + protein + '.csv'`). The code uses `%s` to format the string with the `protein` variable, but it can be written in other ways. 

   - This CSV file is then read into a pandas DataFrame, indicating there is no header, and then writing the header next:

   `mol_ct.columns = ['Number_of_Molecules','file_name']`
   
   - The script calculates how many molecules to sample from each file based on the total number of molecules available by summing the `Number_of_Molecules` column we just created: 
   
   `Total_mols_available = np.sum(mol_ct.Number_of_Molecules)`

   - The number of molecules to sample is calculated as the total number of molecules to sample (`Total_sampling` which =`tot_sampling`) divided by the total number of molecules available (`Total_mols_available`) multiplied by the number of molecules in each file (`elem`): 

   `mol_ct['Sample_for_million'] = [int(Total_sampling/Total_mols_available*elem) for elem in mol_ct.Number_of_Molecules]`

   To me, it is confusing to say that it is `Sample_per_million` when it is really trying to say that this is the proportion of molecules to sample from each file.  
   
   - An updated CSV file (`Mol_ct_file_updated_[project_name].csv`) is saved, including the original counts and the calculated sample sizes per million molecules.
