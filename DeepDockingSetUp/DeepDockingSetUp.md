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

#### `molecular_file_count_updated.py`

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

##### Debugging Print Function
It overrides the built-in `print` function for debugging, prefixing output with "\t molecular_file_count_updated: " to differentiate script output from other messages.

##### Argument Parsing
The five required command-line arguments:
- **project_name** (`-pt`): Name of the DD project.
- **n_iteration** (`-it`): The current iteration number of the DD process.
- **data_directory** (`-cdd`): The path to the directory containing molecule database files.
- **tot_process** (`-t_pos`): The number of CPUs to use for multiprocessing.
- **tot_sampling** (`-t_samp`): The total number of molecules to sample in the current iteration.

All are set from the information we pass in `phase 1.sh`. Importantly, `t_samp` is found via `mols_to_dock`.

##### Main Functionality
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
   
   - An updated CSV file (`Mol_ct_file_updated_[project_name].csv`) is saved, including the original counts and the calculated sample sizes per million molecules (or more plainly the proportion to be sampled from the file in question).

So keep in mind, two files are written when this script is run: 
- `data_directory/Mol_ct_file_[project_name].csv`: The original file with the molecule counts and file names.
- `data_directory/Mol_ct_file_updated_[project_name].csv`: The updated file with the molecule counts, file names, and sample sizes per million molecules.  

And remember as well that `data_directory/` is the `pred_directory` variable we set in `phase_1.sh` and that will update as we move through iterations. 

Now we can get to `sampling.py` which is the next script called in `phase_1.sh`. 

#### `sampling.py`

```sh
python scripts_1/sampling.py -pt $protein -fp $file_path -it $1 -dd $pred_directory -t_pos $t_cpu -tr_sz $mol_to_dock -vl_sz $n_mol
```

```py
from contextlib import closing
from multiprocessing import Pool
import pandas as pd
import numpy as np
import argparse
import glob
import time
import os

try:
    import __builtin__
except ImportError:
    # Python 3
    import builtins as __builtin__

# For debugging purposes only:
def print(*args, **kwargs):
    __builtin__.print('\t sampling: ', end="")
    return __builtin__.print(*args, **kwargs)


parser = argparse.ArgumentParser()
parser.add_argument('-pt', '--project_name',required=True,help='Name of the DD project')
parser.add_argument('-fp', '--file_path',required=True,help='Path to the project directory, excluding project directory name')
parser.add_argument('-it', '--n_iteration',required=True,help='Number of current iteration')
parser.add_argument('-dd', '--data_directory',required=True,help='Path to directory containing the remaining molecules of the database; if first iteration, path to Morgan fingerprints of the database, if other iteration path to morgan_1024_predictions folder of the previous iteration')
parser.add_argument('-t_pos', '--tot_process',required=True,help='Number of CPUs to use for multiprocessing')
parser.add_argument('-tr_sz', '--train_size',required=True,help='Size of training set')
parser.add_argument('-vl_sz', '--val_size',required=True,help='Size of validation and test set')
io_args = parser.parse_args()

protein = io_args.project_name
file_path = io_args.file_path
n_it = int(io_args.n_iteration)
data_directory = io_args.data_directory
tot_process = int(io_args.tot_process)
tr_sz = int(io_args.train_size)
vl_sz = int(io_args.val_size)
rt_sz = tr_sz/vl_sz

print("Parsed Args:")
print(" - Iteration:", n_it)
print(" - Data Directory:", data_directory)
print(" - Training Size:", tr_sz)
print(" - Validation Size:", vl_sz)


def train_valid_test(file_name):
    sampling_start_time = time.time()
    f_name = file_name.split('/')[-1]
    mol_ct = pd.read_csv(data_directory+"/Mol_ct_file_updated_%s.csv"%protein, index_col=1)
    if n_it == 1:
        to_sample = int(mol_ct.loc[f_name].Sample_for_million/(rt_sz+2))
    else:
        to_sample = int(mol_ct.loc[f_name].Sample_for_million/3)

    total_len = int(mol_ct.loc[f_name].Number_of_Molecules)
    shuffle_array = np.linspace(0, total_len-1, total_len)
    seed = np.random.randint(0, 2**32)
    np.random.seed(seed=seed)
    np.random.shuffle(shuffle_array)

    if n_it == 1:
        train_ind = shuffle_array[:int(rt_sz*to_sample)]
        valid_ind = shuffle_array[int(to_sample*rt_sz):int(to_sample*(rt_sz+1))]
        test_ind = shuffle_array[int(to_sample*(rt_sz+1)):int(to_sample*(rt_sz+2))]
    else:
        train_ind = shuffle_array[:to_sample]
        valid_ind = shuffle_array[to_sample:to_sample*2]
        test_ind = shuffle_array[to_sample*2:to_sample*3]

    train_ind_dict = {}
    valid_ind_dict = {}
    test_ind_dict = {}

    train_set = open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/train_set.txt", 'a')
    test_set = open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/test_set.txt", 'a')
    valid_set = open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/valid_set.txt", 'a')

    for i in train_ind:
        train_ind_dict[i] = 1
    for j in valid_ind:
        valid_ind_dict[j] = 1
    for k in test_ind:
        test_ind_dict[k] = 1

    # Opens the file and write the test, train, and valid files
    with open(file_name, 'r') as ref:
        for ind, line in enumerate(ref):
            molecule_id = line.strip().split(',')[0]
            if ind == 1:
                print("molecule_id:", molecule_id)

            # now we write to the train, test, and validation sets
            if ind in train_ind_dict.keys():
                train_set.write(molecule_id + '\n')
            elif ind in valid_ind_dict.keys():
                valid_set.write(molecule_id + '\n')
            elif ind in test_ind_dict.keys():
                test_set.write(molecule_id + '\n')

    train_set.close()
    valid_set.close()
    test_set.close()
    print("Process finished sampling in " + str(time.time()-sampling_start_time))

if __name__ == '__main__':
    try:
        os.mkdir(file_path+'/'+protein+"/iteration_"+str(n_it))
    except OSError:
        pass

    f_names = []
    for f in glob.glob(data_directory+'/*.txt'):
        f_names.append(f)

    t = time.time()
    print("Starting Processes...")
    with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
        pool.map(train_valid_test, f_names)

    print("Compressing smile file...")
    print("Sampling Complete - Total Time Taken:", time.time()-t)
```
##### Debugging Print Function
It overrides the built-in `print` function for debugging, prefixing output with "\t sampling: " to differentiate script output from other messages.

##### Argument Parsing
The seven required command-line arguments:
- **project_name** (`-pt`): Name of the DD project.
- **file_path** (`-fp`): Path to the project directory, excluding project directory name. This is the `file_path` variable we set in `phase_1.sh`.
- **n_iteration** (`-it`): The current iteration number of the DD process. This is the `$1` variable we set in `phase_1.sh`.
- **data_directory** (`-dd`): The path to the directory containing molecule database files. This is the `pred_directory` variable we set in `phase_1.sh`.
- **tot_process** (`-t_pos`): The number of CPUs to use for multiprocessing. This is the `$t_cpu` variable we set in `phase_1.sh`.
- **train_size** (`-tr_sz`): The size of the training set. This is the `$mol_to_dock` variable we set in `phase_1.sh`.
- **val_size** (`-vl_sz`): The size of the validation and test sets. This is the `$n_mol` variable we set in `phase_1.sh`.

Let's remind ourselves about the difference between `$mol_to_dock` and `$n_mol`:

- `$mol_to_dock` (used for *training set* size): comes from the arguments set in the phase_1.sh execution: 

    `mol_to_dock=$5`

    This is the number of molecules to sample from the database for the current iteration. Part of the elegance here is that when you call phase_1.sh, the fifth argument is *just* setting the training size, meaning you can make this any number you want in subsequent iterations (typically in the first they would be equal to what you have in line 8 of the `logs.txt`, but you don't *have* to do that). 

- `$n_mol` (used for *validation and test set* size): comes from the `logs.txt` file:

    `n_mol=`sed -n '8p' $3/$4/logs.txt`
    
    Then the actual amount we will dock will include 2x n_mol (since we need validation and test to be indentical in size) and mol_to_dock if we are in iteration 1, but otherwise we will just be adding to our training set size. See below:

    ```sh
    if [ $1 == 1 ]
    then 
        to_d=$((n_mol+n_mol+mol_to_dock))
    else
        to_d=$mol_to_dock
    fi
    ```

Now, let's return to the python script. After the imports and print debugging function, we get to the main function we will end up using with `multiprocessing`:
```py
def train_valid_test(file_name):
    sampling_start_time = time.time()
    f_name = file_name.split('/')[-1]
    mol_ct = pd.read_csv(data_directory+"/Mol_ct_file_updated_%s.csv"%protein, index_col=1)
    if n_it == 1:
        to_sample = int(mol_ct.loc[f_name].Sample_for_million/(rt_sz+2))
    else:
        to_sample = int(mol_ct.loc[f_name].Sample_for_million/3)
```

This function is called for each file in the `f_names` list, which is created (later in the script) by `glob.glob` to find all `.txt` files in the specified data directory. This is the reason why I often had issues getting files to properly parse, because I'd have `.smi` files in the directory and not `.txt` files.

We get f_name by splitting the file name on the `/` character and taking the last element of the resulting list. This is the file name without the path. This could probably be done more elegantly with `os.path.basename` but this works. 

Next, we read in the `Mol_ct_file_updated_[project_name].csv` file we created in `molecular_file_count_updated.py` and set the index to the file name. This is important because we will be using the file name to look up the number of molecules to sample from the file via the `loc` method (`mol_ct.loc[f_name].Sample_for_million`). I have not seen this before, but this is a clever way to get the value we want, which is the "sample proportion" (called Sample_for_million) from the file we are currently processing. 

If we are in iteration 1, we set `to_sample` as the sample proportion dividided `rt_sz`, which is defined as `rt_sz = tr_sz/vl_sz`.

```py
    total_len = int(mol_ct.loc[f_name].Number_of_Molecules)
    shuffle_array = np.linspace(0, total_len-1, total_len)
    seed = np.random.randint(0, 2**32)
    np.random.seed(seed=seed)
    np.random.shuffle(shuffle_array)
```
I think it's helpful to see what these files look like in an example before I continue explaining the code to myself. Let's use my 3_EnamineDiverse project as an example. In `interation_1/morgan_1024_predictions/` I have: 

```sh
Mol_ct_file_3_EnamineDiverse.csv
Mol_ct_file_updated_3_EnamineDiverse.csv
part_00_smiles.txt
...
part_27_smiles.txt
passed_file_ct.txt
```
We are looking at: `Mol_ct_file_updated_3_EnamineDiverse.csv` right now. 

It looks like this: 

```py
Number_of_Molecules,file_name,Sample_for_million
14143753,part_10_smiles.txt,213671
```
NOTE: You write 'back in time' here, so Iteration_1/morgan_1024_predictions/Mol_ct_updated_3_EnamineDiverse.csv is writen during phase_1.sh in *iteration 2*. So when I discuss the above file as example, remember the sampling is going to go into iteration 2, not iteration 1. Iteration 1 instead writes to the library data directory. This is most obvious when you `fd` for the file: 

```sh
fd Mol_ct_file_updated_3_EnamineDiverse.csv
enamine/rename/enumerate/cat_for_isomer_deletion/enamine_real_28_fp/Mol_ct_file_updated_3_EnamineDiverse.csv
```

Anyway, the `Mol_ct_file_updated_3_EnamineDiverse.csv` file is created in `molecular_file_count_updated.py` and is read in here. The `f_name` variable is the file name without the path. We use the number of molecules to construct the array that we will then shuffle with a random seed. Relatively straight forward if you get where the paths are going. 

If you're in the first iteration, the bounds of the array are taken from the shuffled array and plit by the `rt_sz` variable. I am still a little unsure of why it is set up this way, but the operation itself makes sense. 

If you're not in the first iteration, then the bounds of the array are taken from the shuffled array and split into three equal parts. 

Regardless of iteration however, you always get a train/valid/test split. It just that in later iterations, all the sets will contribute to making the training set.

Then, the *_set.txt is saved in the iteration_X directory and written.

I find the part with the indices setting hard to follow, so here is one potential answer:

`shuffle_array`: This is a numpy array that starts as a sequence of integers from 0 to total_len-1. After shuffling, the order of these integers is randomized. For example, if `total_len` is 5, `shuffle_array` might start as `[0, 1, 2, 3, 4]` and then become something like `[3, 0, 4, 1, 2]` after shuffling.

`train_ind`, `valid_ind`, `test_ind`: These are slices of the shuffled array, so they are also numpy arrays containing a subset of the integers from 0 to total_len-1, in a random order. For example, if `to_sample` is 2, `train_ind` might be `[3, 0]`, `valid_ind` might be `[4, 1]`, and `test_ind` might be [2].

`train_ind_dict`, `valid_ind_dict`, `test_ind_dict`: These are dictionaries where the keys are the integers in `train_ind`, `valid_ind`, and `test_ind`, respectively, and the values are all 1. For example, if `train_ind` is `[3, 0]`, `train_ind_dict` would be `{3: 1, 0: 1}`.

Then, when we act on the file_name (the *.txt files in the prediction directory) we use the index of the file to sample. So using `train_ind_dict` = `{3: 1, 0: 1}`

and a file containing the following, the bold would be written to the train_set.txt file:

[0] **Z1449657345_Isomer1,0.998218834400177**  
[1] Z2861190319_Isomer1,0.04097618907690048  
[2] Z2861190319_Isomer2,0.2449219524860382  
[3] **Z6127089457_Isomer1,0.24301858246326447**  
[4] Z6127400168_Isomer1,0.7136911153793335  

There's a couple of things here that are a bit confusing. One is that while I originally thought having random seeds per process would cause issues, it doesnt if the `*.txt` files are unique (then any sort of indice issue doesnt matter, you can accidentally write the same line to two different datasets if every line is unique).

The other is that it doesn't seem like you actually need a dictionary. It seems like a dictionary here is arbitrarily set to have key, value = indice, 1; but `1` is constant. This is likely because dictionary membership is faster to check than list membership (for whatever reason), one suggestion I found was to use set() instead, like so: 

```py
train_ind_set = set(train_ind)
valid_ind_set = set(valid_ind)
test_ind_set = set(test_ind)

# ...

if ind in train_ind_set:
    train_set.write(molecule_id + '\n')
elif ind in valid_ind_set:
    valid_set.write(molecule_id + '\n')
elif ind in test_ind_set:
    test_set.write(molecule_id + '\n')
```

I'm not sure if this is faster, but it's a good idea to keep in mind. 

However, now we know this script will sucessfully write the `test_set.txt`, `train_set.txt`, and `valid_set.txt` files. 

Next we will look at the `sanity_check.py` script, which ensures deduplication.

```sh
python scripts_1/sanity_check.py -pt $protein -fp $file_path -it $1
```
```py
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-pt','--project_name',required=True,help='Name of project')
parser.add_argument('-fp','--file_path',required=True,help='Path to project folder without name of project folder')
parser.add_argument('-it','--n_iteration',required=True,help='Number of current iteration')

io_args = parser.parse_args()
import time

protein = io_args.project_name
file_path = io_args.file_path
n_it = int(io_args.n_iteration)

old_dict = {}
for i in range(1,n_it):
    with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(i)+'/training_labels*')[-1]) as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split(',')[-1]
            old_dict[tmpp] = 1
    with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(i)+'/validation_labels*')[-1]) as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split(',')[-1]
            old_dict[tmpp] = 1
    with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(i)+'/testing_labels*')[-1]) as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split(',')[-1]
            old_dict[tmpp] = 1

t=time.time()
new_train = {}
new_valid = {}
new_test = {}
with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(n_it)+'/train_set*')[-1]) as ref:
    for line in ref:
        tmpp = line.strip().split(',')[0]
        new_train[tmpp] = 1
with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(n_it)+'/valid_set*')[-1]) as ref:
    for line in ref:
        tmpp = line.strip().split(',')[0]
        new_valid[tmpp] = 1
with open(glob.glob(file_path+'/'+protein+'/iteration_'+str(n_it)+'/test_set*')[-1]) as ref:
    for line in ref:
        tmpp = line.strip().split(',')[0]
        new_test[tmpp] = 1
print(time.time()-t)

t=time.time()
for keys in new_train.keys():
    if keys in new_valid.keys():
        new_valid.pop(keys)
    if keys in new_test.keys():
        new_test.pop(keys)
for keys in new_valid.keys():
    if keys in new_test.keys():
        new_test.pop(keys)
print(time.time()-t)

for keys in old_dict.keys():
    if keys in new_train.keys():
        new_train.pop(keys)
    if keys in new_valid.keys():
        new_valid.pop(keys)
    if keys in new_test.keys():
        new_test.pop(keys)
        
with open(file_path+'/'+protein+'/iteration_'+str(n_it)+'/train_set.txt','w') as ref:
    for keys in new_train.keys():
        ref.write(keys+'\n')
with open(file_path+'/'+protein+'/iteration_'+str(n_it)+'/valid_set.txt','w') as ref:
    for keys in new_valid.keys():
        ref.write(keys+'\n')
with open(file_path+'/'+protein+'/iteration_'+str(n_it)+'/test_set.txt','w') as ref:
    for keys in new_test.keys():
        ref.write(keys+'\n')
```

This script is pretty straight forward. It takes the project name, file path, and iteration number as arguments. Then it creates a dictionary of all the molecules that have been used in previous iterations via the `old_dict` dictionary up to the current iteration via `for i in range(1,n_it)`. 

The extraction logic uses: 

`tmpp = line.strip().split(',')[-1]`

To set the keys of the dictionary to the molecule ID. Notably, we are checking the `*_labels.txt` files, not the `*_set.txt` files when looking to previous iterations.

This is the structure of the `*_labels.txt` files: 

```
r_i_docking_score,ZINC_ID
-10.4494,PV-006642569130_Isomer1
```
As can be seen, the last element (found via `[-1]`), is responsible for finding the molecule ID.

The "new" keys are then from the `*_set.txt` files. 

We can then use the `pop` method to remove any molecules that are in the `new` dictionaries from the `old` dictionary. This is done for the `train`, `valid`, and `test` dictionaries. 

Finally, the `*_set.txt` files are written with the remaining molecules, ensuring previous labels have not been re-used. 

We can now move on to `extracting_morgan.py`:

#### `extracting_morgan.py`

```sh
python scripts_1/extracting_morgan.py -pt $protein -fp $file_path -it $1 -md $morgan_directory -t_pos $t_cpu
```

```py
# Reads the ids found in sampling and finds the corresponding morgan fingerprint
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument('-pt', '--project_name', required=True, help='Name of the DD project')
parser.add_argument('-fp', '--file_path', required=True, help='Path to the project directory, excluding project directory name')
parser.add_argument('-it', '--n_iteration', required=True, help='Number of current iteration')
parser.add_argument('-md', '--morgan_directory', required=True, help='Path to directory containing Morgan fingerprints for the database')
parser.add_argument('-t_pos', '--tot_process', required=True, help='Number of CPUs to use for multiprocessing')

io_args = parser.parse_args()

import os
from multiprocessing import Pool
import time
from contextlib import closing
import numpy as np

protein = io_args.project_name
file_path = io_args.file_path
n_it = int(io_args.n_iteration)
morgan_directory = io_args.morgan_directory
tot_process = int(io_args.tot_process)


def extract_morgan(file_name):
    train = {}
    test = {}
    valid = {}
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/train_set.txt", 'r') as ref:
        for line in ref:
            train[line.rstrip()] = 0
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/valid_set.txt", 'r') as ref:
        for line in ref:
            valid[line.rstrip()] = 0
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/test_set.txt", 'r') as ref:
        for line in ref:
            test[line.rstrip()] = 0

    # for file_name in file_names:
    ref1 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'train_' + file_name.split('/')[-1], 'w')
    ref2 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'valid_' + file_name.split('/')[-1], 'w')
    ref3 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'test_' + file_name.split('/')[-1], 'w')

    with open(file_name, 'r') as ref:
        flag = 0
        for line in ref:
            tmpp = line.strip().split(',')[0]
            if tmpp in train.keys():
                train[tmpp] += 1
                fn = 1
                if train[tmpp] == 1: flag = 1
            elif tmpp in valid.keys():
                valid[tmpp] += 1
                fn = 2
                if valid[tmpp] == 1: flag = 1
            elif tmpp in test.keys():
                test[tmpp] += 1
                fn = 3
                if test[tmpp] == 1: flag = 1
            if flag == 1:
                if fn == 1:
                    ref1.write(line)
                if fn == 2:
                    ref2.write(line)
                if fn == 3:
                    ref3.write(line)
            flag = 0


def alternate_concat(files):
    to_return = []
    with open(files, 'r') as ref:
        for line in ref:
            to_return.append(line)
    return to_return


def delete_all(files):
    os.remove(files)


def morgan_duplicacy(f_name):
    flag = 0
    mol_list = {}
    ref1 = open(f_name[:-4] + '_updated.csv', 'a')
    with open(f_name, 'r') as ref:
        for line in ref:
            tmpp = line.strip().split(',')[0]
            if tmpp not in mol_list:
                mol_list[tmpp] = 1
                flag = 1
            if flag == 1:
                ref1.write(line)
                flag = 0
    os.remove(f_name)


if __name__ == '__main__':
    try:
        os.mkdir(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan')
    except:
        pass

    files = []
    for f in glob.glob(morgan_directory + "/*.txt"):
        files.append(f)

    t = time.time()
    with closing(Pool(np.min([tot_process, len(files)]))) as pool:
        pool.map(extract_morgan, files)
    print(time.time() - t)

    all_to_delete = []
    for type_to in ['train', 'valid', 'test']:
        t = time.time()
        files = []
        for f in glob.glob(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + type_to + '*'):
            files.append(f)
            all_to_delete.append(f)
        print(len(files))
        if len(files) == 0:
            print("Error in address above")
            break
        with closing(Pool(np.min([tot_process, len(files)]))) as pool:
            to_print = pool.map(alternate_concat, files)
        with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + type_to + '_morgan_1024.csv',
                  'w') as ref:
            for file_data in to_print:
                for line in file_data:
                    ref.write(line)
        to_print = []
        print(type_to, time.time() - t)

    f_names = []
    for f in glob.glob(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/*morgan*'):
        f_names.append(f)

    t = time.time()
    with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
        pool.map(morgan_duplicacy, f_names)
    print(time.time() - t)

    with closing(Pool(np.min([tot_process, len(all_to_delete)]))) as pool:
        pool.map(delete_all, all_to_delete)
```

##### Functions: 

We are familiar with imports now, so let's move to function-by-function explanations: 

```py
def extract_morgan(file_name):
    train = {}
    test = {}
    valid = {}
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/train_set.txt", 'r') as ref:
        for line in ref:
            train[line.rstrip()] = 0
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/valid_set.txt", 'r') as ref:
        for line in ref:
            valid[line.rstrip()] = 0
    with open(file_path + '/' + protein + "/iteration_" + str(n_it) + "/test_set.txt", 'r') as ref:
        for line in ref:
            test[line.rstrip()] = 0
```
Here we parse the `*_set.txt` files into dictionaries. The `rstrip()` method removes the trailing newline character. 

So, if we consider this as the structure of some `*_set.txt` file: 

PV-004170659562_Isomer1
PV-004170659562_Isomer2

The dictionary would look like this: 

```py
{
    "PV-004170659562_Isomer1": 0,
    "PV-004170659562_Isomer2": 0
}
```

The reason for the `rstrip()` method is that otherwise we would include the trailing newline character (`\n` in python).


``` py
    # for file_name in file_names:
    ref1 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'train_' + file_name.split('/')[-1], 'w')
    ref2 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'valid_' + file_name.split('/')[-1], 'w')
    ref3 = open(
        file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'test_' + file_name.split('/')[-1], 'w')
```

Here we are using the `file_name` argument that is passed to the `extract_morgan(file_name)` function. In the `main()` function, we call it like so: 

```py
with closing(Pool(np.min([tot_process, len(files)]))) as pool:
    pool.map(extract_morgan, files)
```

So, `file_name` is a file path to a file that contains the morgan fingerprints. They are found like so, also in `main()`:

```py
files = []
for f in glob.glob(morgan_directory + "/*.txt"):
    files.append(f)
```

The `morgan_directory` variable comes via `argparse` and is defined: 

```py
parser.add_argument('-md', '--morgan_directory', required=True, help='Path to directory containing Morgan fingerprints for the database')
```
So, `morgan_directory` is the path to the directory that contains the morgan fingerprints, which end in `.txt`. It seems like we are going to spawn many `train_`, `valid_`, and `test_` files, as `file_name` is a list of file paths that would look like this: 

```sh
path/to/morgan_directory/part_00_smiles.txt
...
path/to/morgan_directory/part_27_smiles.txt
```
So, the `file_name.split('/')[-1]` would look like this: 

```sh
part_00_smiles.txt
...
part_27_smiles.txt
```
But, keep in mind, these files are being opened in:

`ref1 = open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + 'train_' + file_name.split('/')[-1], 'w')`

So this is in the `morgan` directory of the current iteration.

Now the final part of the function:

```py
    with open(file_name, 'r') as ref:
        flag = 0
        for line in ref:
            tmpp = line.strip().split(',')[0]
            if tmpp in train.keys():
                train[tmpp] += 1
                fn = 1
                if train[tmpp] == 1: flag = 1
            elif tmpp in valid.keys():
                valid[tmpp] += 1
                fn = 2
                if valid[tmpp] == 1: flag = 1
            elif tmpp in test.keys():
                test[tmpp] += 1
                fn = 3
                if test[tmpp] == 1: flag = 1
            if flag == 1:
                if fn == 1:
                    ref1.write(line)
                if fn == 2:
                    ref2.write(line)
                if fn == 3:
                    ref3.write(line)
            flag = 0
```
We open the `file_name` we've discussed and iterate over it. We are looking for the molecule ID, which is the first element of the line (because it is split by `,`). The `.strip()` is for leading/trailing whitespace. We set a flag to 0 and then check if the molecule ID is in the `train`, `valid`, or `test` dictionaries. If it is, we increment the value of the dictionary by 1. We also set `fn` to 1, 2, or 3 depending on which dictionary it is in. If the value of the dictionary is 1, we set the flag to 1. 

The trick to understanding this is to realize that the entire set of if statements has to run. So, if the molecule ID is in the `train` dictionary, we increment the value by 1, set `fn` to 1, and then check if the value is 1. If it is, we set the flag to 1. Then, we get to if `flag == 1`, which it is if we have a match, then we also know via `fn` which file to write to. Then we reset `flag` to 0 and go to the next line. If it's not matched, then flag is never changed from 0, and therefor, nothing is ever written. If the watch is in a different dictionary, then we just write to the `fn` number it corresponds to. 

The output, presumably, are going to be files that looks like this: 

```sh
path/to/morgan/directory/train_part_00_smiles.txt
```

And so on for `valid_`, `test_` prefixes and for the `file_names` we are acting on in `main()`.

We are writing the entire line of `morgan_library_fp.txt` file we are reading to the file, so the output would look like this, using an example I have: 

```sh
head -n 1 part_00_smiles.txt
```

```sh
Z6430701089_Isomer1,3,15,24,33,36,59,64,80,90,123,128,188,197,222,270,282,284,356,367,371,375,389,394,422,441,447,546,547,606,614,650,656,658,689,698,703,726,807,849,854,893,899,926,975,1019,1021
```

Now, the next steps are likely to appropriately concat these files and remove the the many files we just made that we won't need. There also appears to be a step to check for duplications, though I am not sure why this is necessary. 

Here are the last three functions before we begin using them in `main()`:

```py
def alternate_concat(files):
    to_return = []
    with open(files, 'r') as ref:
        for line in ref:
            to_return.append(line)
    return to_return


def delete_all(files):
    os.remove(files)


def morgan_duplicacy(f_name):
    flag = 0
    mol_list = {}
    ref1 = open(f_name[:-4] + '_updated.csv', 'a')
    with open(f_name, 'r') as ref:
        for line in ref:
            tmpp = line.strip().split(',')[0]
            if tmpp not in mol_list:
                mol_list[tmpp] = 1
                flag = 1
            if flag == 1:
                ref1.write(line)
                flag = 0
    os.remove(f_name)
```

The `alternate_concat` function is used to read the files we just made. It is used in `main()` like so: 

```py
with closing(Pool(np.min([tot_process, len(files)]))) as pool:
    to_print = pool.map(alternate_concat, files)
```

Though I don't really get why? As we seem to overwrite to_print in the next step. 

The `delete_all` function is used to remove the files we just made. It is used in `main()` like so: 

```py
with closing(Pool(np.min([tot_process, len(all_to_delete)]))) as pool:
    pool.map(delete_all, all_to_delete)
```
And finally, the `morgan_duplicacy` function is used to check for duplications. It is used in `main()` like so: 

```py
with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
    pool.map(morgan_duplicacy, f_names)
```
We start with `flag = 0` and an empty dictionary `mol_list`. We are going to strip the file extension from `f_name` via `f_name[:-4]`, add `'_updated.csv'` to it, and open it for appending. Then we open `f_name` and iterate over it. We are looking for the molecule ID, which is the first element of the line (because it is split by `,`). The `.strip()` is for leading/trailing whitespace. We set a flag to 0 and then check if the molecule ID is in the `mol_list` dictionary. If it is not, we add it to the dictionary and set the flag to 1. If the flag is 1, we write the line to the file we opened. Then we reset the flag to 0. This is pretty similar to what we did before. 

The part that is confusing me is what f_name is, because if I expected it to be something like `train_part_00_smiles.txt`, but it is not. It would have to be something like `train.txt` because our output files in `/mnt/data/dk/work/DeepDocking/projects/3_EnamineDiverse/iteration_1/morgan` are: 

```sh
test_morgan_1024_updated.csv
train_morgan_1024_updated.csv
valid_morgan_1024_updated.csv
```
Ah, so what's happening is further down in `main()`:

```py
    all_to_delete = []
    for type_to in ['train', 'valid', 'test']:
        t = time.time()
        files = []
        for f in glob.glob(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + type_to + '*'):
            files.append(f)
            all_to_delete.append(f)
        print(len(files))
        if len(files) == 0:
            print("Error in address above")
            break
        with closing(Pool(np.min([tot_process, len(files)]))) as pool:
            to_print = pool.map(alternate_concat, files)
        with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/' + type_to + '_morgan_1024.csv',
                  'w') as ref:
            for file_data in to_print:
                for line in file_data:
                    ref.write(line)
        to_print = []
        print(type_to, time.time() - t)
```
The first part of the code is going to find the `train_`, `test_`, and `valid_` prefixed files we made earlier with `extract_morgan()` using the `glob.glob()` search. 

Then we use `alternate_concat()`. It reads the files we just made (`files = []`) and returns them as a list. So, `to_print` is a list of lists. The first list is the contents of `train_part_00_smiles.txt`, the second list is the contents of `train_part_01_smiles.txt`, and so on. This is why we need this double for loop to write the lines to the file:

```py
for file_data in to_print:
    for line in file_data:
        ref.write(line)
```
This is another detail about python's multiprocessing: 

"The `multiprocessing.Pool.map()` function always returns a list. It applies the function you give it to each element in the iterable you pass, and collects the results in a list. The order of the results in the list corresponds to the order of the elements in the input iterable.

So, if alternate_concat returns a list, `pool.map(alternate_concat, files)` will return a list of lists. If `alternate_concat` returned a dictionary, `pool.map(alternate_concat, files)` would return a list of dictionaries.

In other words, `Pool.map()` doesn't change the type of the return value of the function it's applying. It just collects those return values in a list. The type of the elements in that list depends on what the function returns."

Now, I believe we have solved the question from earlier.

```py
f_names = []
    for f in glob.glob(file_path + '/' + protein + '/iteration_' + str(n_it) + '/morgan/*morgan*'):
        f_names.append(f)

    t = time.time()
    with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
        pool.map(morgan_duplicacy, f_names)
    print(time.time() - t)

    with closing(Pool(np.min([tot_process, len(all_to_delete)]))) as pool:
        pool.map(delete_all, all_to_delete)
```

f_names here is going to be a list of files that look like this: 

```sh
train_morgan_1024.csv
valid_morgan_1024.csv
test_morgan_1024.csv
```
When we run `morgan_duplicacy()` on these files, we are going to check for duplications in the files we just made, remove them, and append `_updated.csv` to the end of the file name.

We are also going to store the `X_morgan_1024.csv` files in a list and pass them to a `delete_all` function, cleaning up the directory. 

So, the final output of this script is going to be: 

```sh
train_morgan_1024_updated.csv
valid_morgan_1024_updated.csv
test_morgan_1024_updated.csv
```

Now, we need to extract the smiles to finish phase_1!

#### `extracting_smiles.py`

```sh
python scripts_1/extracting_smiles.py -pt $protein -fp $file_path -it $1 -smd $smile_directory -t_pos $t_cpu
```

```py
from multiprocessing import Pool
from functools import partial
from contextlib import closing
import argparse
import numpy as np
import pickle
import glob
import time
import os

parser = argparse.ArgumentParser()
parser.add_argument('-pt', '--project_name', required=True, help='Name of the DD project')
parser.add_argument('-fp', '--file_path', required=True, help='Path to the project directory, excluding project directory name')
parser.add_argument('-it', '--n_iteration', required=True, help='Number of current iteration')
parser.add_argument('-smd', '--smile_directory', required=True, help='Path to SMILES directory of the database')
parser.add_argument('-t_pos', '--tot_process', required=True, help='Number of CPUs to use for multiprocessing')

io_args = parser.parse_args()
protein = io_args.project_name
file_path = io_args.file_path
n_it = int(io_args.n_iteration)
smile_directory = io_args.smile_directory
tot_process = int(io_args.tot_process)

# This path is used often so we declare it as a constant here.
ITER_PATH = file_path + '/' + protein + '/iteration_' + str(n_it)

def extract_smile(file_name, train, valid, test):
    # This function extracts the smiles from a file to write them to train, test, and valid files for model training.
    ref1 = open(ITER_PATH + '/smile/' + 'train_' + file_name.split('/')[-1], 'w')
    ref2 = open(ITER_PATH + '/smile/' + 'valid_' + file_name.split('/')[-1], 'w')
    ref3 = open(ITER_PATH + '/smile/' + 'test_' + file_name.split('/')[-1], 'w')

    with open(file_name, 'r') as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split()[-1]
            if tmpp in train.keys():
                train[tmpp] += 1
                if train[tmpp] == 1: ref1.write(line)

            elif tmpp in valid.keys():
                valid[tmpp] += 1
                if valid[tmpp] == 1: ref2.write(line)

            elif tmpp in test.keys():
                test[tmpp] += 1
                if test[tmpp] == 1: ref3.write(line)

def alternate_concat(files):
    # Returns a list of the lines in a file
    with open(files, 'r') as ref:
        return ref.readlines()

def smile_duplicacy(f_name):
    # removes duplicate molec from the file 
    mol_list = {} # keeping track of which mol have been written
    ref1 = open(f_name[:-4] + '_updated.smi', 'a')
    with open(f_name, 'r') as ref:
        for line in ref:
            tmpp = line.strip().split()[-1]
            if tmpp not in mol_list: # avoiding duplicates
                mol_list[tmpp] = 1 
                ref1.write(line)
    os.remove(f_name)

def delete_all(files):
    os.remove(files)

if __name__ == '__main__':
    try:
        os.mkdir(ITER_PATH + '/smile')
    except: # catching exception for when the folder already exists
        pass

    files_smiles = [] # Getting the path to every smile file from docking
    for f in glob.glob(smile_directory + "/*.txt"):
        files_smiles.append(f)

    # return_mols_per_file = []
    # for j in range(3):
    #     ct = 0
    #     for i in range(len(return_mols_per_file)):
    #         ct += len(return_mols_per_file[i][j])
    #     print(ct)

    # for k in range(3):
    #     t = time.time()
    #     for i in range(len(return_mols_per_file)):
    #         for j in range(i + 1, len(return_mols_per_file)):
    #             for keys in return_mols_per_file[i][k].keys():
    #                 if keys in return_mols_per_file[j][k]:
    #                     return_mols_per_file[j][k].pop(keys)
    #     print(time.time() - t)

    # for j in range(3):
    #     ct = 0
    #     for i in range(len(return_mols_per_file)):
    #         ct += len(return_mols_per_file[i][j])
    #     print(ct)

    # train = {}
    # valid = {}
    # test = {}
    # for j in range(3):
    #     for i in range(len(return_mols_per_file)):
    #         for keys in return_mols_per_file[i][j]:
    #             if j == 0:
    #                 train[keys] = 0
    #             elif j == 1:
    #                 valid[keys] = 0
    #             elif j == 2:
    #                 test[keys] = 0

    all_train = {}
    all_valid = {}
    all_test = {}
    with open(ITER_PATH + "/train_set.txt", 'r') as ref:
        for line in ref:
            all_train[line.rstrip()] = 0

    with open(ITER_PATH + "/valid_set.txt", 'r') as ref:
        for line in ref:
            all_valid[line.rstrip()] = 0

    with open(ITER_PATH + "/test_set.txt", 'r') as ref:
        for line in ref:
            all_test[line.rstrip()] = 0
    
    # ? I don't think this does anything either,
    # all these dictionaires are empty when initalized
    # for keys in train.keys():
    #     all_train.pop(keys)

    # for keys in valid.keys():
    #     all_valid.pop(keys)

    # for keys in test.keys():
    #     all_test.pop(keys)

    print(len(all_train), len(all_valid), len(all_test))

    t = time.time()
    with closing(Pool(np.min([tot_process, len(files_smiles)]))) as pool:
        pool.map(partial(extract_smile, train=all_train, valid=all_valid, test=all_test), files_smiles)
    print(time.time() - t)

    all_to_delete = []
    for type_to in ['train', 'valid', 'test']:
        t = time.time()
        files = []
        for f in glob.glob(ITER_PATH + '/smile/' + type_to + '*'):
            files.append(f)
            all_to_delete.append(f)
        print(len(files))
        if len(files) == 0:
            print("Error in address above")
            break
        with closing(Pool(np.min([tot_process, len(files)]))) as pool:
            to_print = pool.map(alternate_concat, files)
        with open(ITER_PATH + '/smile/' + type_to + '_smiles_final.csv', 'w') as ref:
            for file_data in to_print:
                for line in file_data:
                    ref.write(line)
        to_print = []
        print(type_to, time.time() - t)

    f_names = []
    for f in glob.glob(ITER_PATH + '/smile/*final*'):
        f_names.append(f)

    t = time.time()
    with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
        pool.map(smile_duplicacy, f_names)
    print(time.time() - t)

    with closing(Pool(np.min([tot_process, len(all_to_delete)]))) as pool:
        pool.map(delete_all, all_to_delete)
```

##### Extra Code?

There is a code section I commented out because it doesn't seem to do anything. I am not sure why it is there. 

##### Functions

These functions are very similar to `extracting_morgan.py`, which makes lot of sense. At first glance, it is written a little more cleanly, for example, using `ITER_PATH` as a constant for the path to the iteration directory. There are other changes as well, again, they seem to be synonmous functionally but I'll investigate to be more sure.

The `extract_smile()` function is used to extract the smiles from a file to write them to train, test, and valid files for model training. It has some key differences from `extract_morgan()` however, most notably in that it is taking more arguments (perhaps unnessecary complexity).

```py
def extract_smile(file_name, train, valid, test):
```
Which is different from `extract_morgan()`

```py
def extract_morgan(file_name):
```

The `extract_smile()` function is used in `main()` like so: 

```py
with closing(Pool(np.min([tot_process, len(files_smiles)]))) as pool:
    pool.map(partial(extract_smile, train=all_train, valid=all_valid, test=all_test), files_smiles)
```
The really notable difference here is the use of the `partial()` function, from `functools`. Here is an explanation I found:

>"The `extract_smile` function takes four arguments: `file_name`, `train`, `valid`, and `test`. It opens a file with the given file_name, reads each line, and checks if the last element of the line (after splitting by spaces) is in the train, valid, or test dictionaries. If it is, it increments the corresponding count in the dictionary and writes the line to the corresponding file if it's the first occurrence.
>
>When you use `partial(extract_smile, train=train, valid=valid, test=test)`, you're creating a new function that has the train, valid, and test arguments pre-filled. This new function only needs one argument: `file_name`.
>
>This is useful when you use `pool.map` to apply the function to each element in `files_smiles`. `pool.map` only provides one argument to the function (each element in `files_smiles`), so by using `partial`, you can create a version of `extract_smile` that only needs one argument. This allows you to use `extract_smile` with `pool.map` without having to modify the function or the way you call `pool.map`.
>
>In other words, partial allows you to "adapt" the `extract_smile` function to be used with `pool.map` in a multiprocessing context."

Let's continue look at the `extract_smile()` function. 

```py
def extract_smile(file_name, train, valid, test):
    # This function extracts the smiles from a file to write them to train, test, and valid files for model training.
    ref1 = open(ITER_PATH + '/smile/' + 'train_' + file_name.split('/')[-1], 'w')
    ref2 = open(ITER_PATH + '/smile/' + 'valid_' + file_name.split('/')[-1], 'w')
    ref3 = open(ITER_PATH + '/smile/' + 'test_' + file_name.split('/')[-1], 'w')

    with open(file_name, 'r') as ref:
        ref.readline()
        for line in ref:
            tmpp = line.strip().split()[-1]
            if tmpp in train.keys():
                train[tmpp] += 1
                if train[tmpp] == 1: ref1.write(line)

            elif tmpp in valid.keys():
                valid[tmpp] += 1
                if valid[tmpp] == 1: ref2.write(line)

            elif tmpp in test.keys():
                test[tmpp] += 1
                if test[tmpp] == 1: ref3.write(line)
```

You can see that although the `train`, `valid`, and `test` arguments are used - they are not initalized here. They require the previous portion of `main()` that sets them as empty dictionaries:

```py
all_train = {}
all_valid = {}
all_test = {}
```
And then populates them with the `train_set.txt`, `valid_set.txt`, and `test_set.txt` files.

```py
with open(ITER_PATH + "/train_set.txt", 'r') as ref:
        for line in ref:
            all_train[line.rstrip()] = 0

    with open(ITER_PATH + "/valid_set.txt", 'r') as ref:
        for line in ref:
            all_valid[line.rstrip()] = 0

    with open(ITER_PATH + "/test_set.txt", 'r') as ref:
        for line in ref:
            all_test[line.rstrip()] = 0
```

There is then logic to remove the keys from the dictionaries that are already in the `train`, `valid`, and `test` sets, but those are empty dictionaries so it doesn't do anything. 

Returning to the function, we skip the first line with readline() and then iterate through the file (`file_name`) that is called on the library of SMILES (`files_smiles`). Note: I don't think this `readline()` approach is a great idea, at least without warning the user. Looking at my SMILES library, my part_00_smiles.txt does not have a header. The presence of this header, or the requirement for it, is very software dependent (e.g. does MayaChemTools or Schrodinger require the header? When they output files, do they include the header?). If you tell the user that their SMILES files (which are also called .txt files, which I don't love) need a header, then they can add it and this script can remain mostly as is. 

Here is my `part_00_smiles.txt` file's first line:

Cc1cccc2c1OC[C@@H]2C(=O)NCC1(CO)CNC1 Z6430701089_Isomer1

We can see the last element is the ID. Using `tmpp = line.strip().split()[-1]` is splitting on ` ` (space) and taking the last element. This is a little dangerous, because if the ID has a space in it, it will be split into two elements. I don't think this is a problem for the current data, but it is something to be aware of. 

The `if` statements are checking if the ID is in the `train`, `valid`, or `test` dictionaries. If it is, it increments the count in the dictionary and writes the line to the corresponding file if it's the first occurrence. Since the the dictionaries were previously populated and all if statements are sequential, this should only write the first occurrence of each ID to the corresponding file. 

The `alternate_concat()` function is used to return a list of the lines in a file. It is used in `main()` like so: 
```py
with closing(Pool(np.min([tot_process, len(files)]))) as pool:
    to_print = pool.map(alternate_concat, files)
```

And it is defined like so: 
```py
def alternate_concat(files):
    # Returns a list of the lines in a file
    with open(files, 'r') as ref:
        return ref.readlines()
```

This is functionally equivalent to the `extract_morgan()` function's `alternate_concat()` function. Though, I think the `extract_morgan()` function's `alternate_concat()` function is a little more clear and potentially more memory safe. 

The `smile_duplicacy()` function is used to remove duplicate molecules from the file. It is used in `main()` like so: 
```py
with closing(Pool(np.min([tot_process, len(f_names)]))) as pool:
    pool.map(smile_duplicacy, f_names)
```
As far as I can tell, it behaves just like `morgan_duplicacy()` from `extracting_morgan.py`. 

The `delete_all()` function is used to delete all files. It is used in `main()` like so: 
```py
with closing(Pool(np.min([tot_process, len(all_to_delete)]))) as pool:
    pool.map(delete_all, all_to_delete)
```

Also, as far as I can tell, it behaves just like `delete_all()` from `extracting_morgan.py`. 

The `main()` function is very similar to `extracting_morgan.py`'s `main()` function, except for these likely extraneous code blocks and the use of `partial()` in `pool.map()`. 

Now, we have the SMILES files in the `ITER_PATH + '/smile/'` directory. 

### Phase 2 and Phase 3

Both of these phases are relatively straight forward and control the ligand preparation and docking you do. I'll skip it here for brevity.

It does give me the idea of hosting a version of Deep Docking on my Github that is configured towards my group's use (e.g. LigPrep/Glide).

I could also store the environment necessary to interact with our GPUs. Maybe update to Poetry? 

### Phase 4

Phase 4 extracts the labels we've generated in Phase 2/3 and puts them into a format that can be used for training. 

It is run like so: 

```sh
sbatch phase_4.sh current_iteration 3 path_project project name_gpu_partition tot_number_iterations percent_first_mols percent_last_mols_value recall_value 00-15:00 conda_env
```

Where `current_iteration`  ($1) is the current iteration, `3` ($2) is the total amount of processors available, `path_project` ($3) is the path to the project, `project` ($4) is the name of the project, `name_gpu_partition` ($5) is the name of the partition, `tot_number_iterations` ($6) is the total number of iterations, `percent_first_mols` ($7) is the percentage of the first molecules to be used, `percent_last_mols_value` ($8) is the percentage of the last molecules to be used, `recall_value` ($9) is the recall value, `00-15:00` ($10) is the time, and `conda_env` ($11) is the conda environment. 


```sh
#!/bin/bash
#SBATCH --cpus-per-task=3
#SBATCH --ntasks=1
#SBATCH --mem=0               # memory per node
#SBATCH --job-name=phase_4

env=${11}
time=${10}

source ~/.bashrc
conda activate $env

file_path=`sed -n '1p' $3/$4/logs.txt`
protein=`sed -n '2p' $3/$4/logs.txt`

morgan_directory=`sed -n '4p' $3/$4/logs.txt`
smile_directory=`sed -n '5p' $3/$4/logs.txt`
nhp=`sed -n '7p' $3/$4/logs.txt`    # number of hyperparameters
sof=`sed -n '6p' $3/$4/logs.txt`    # The docking software used

rec=$9

num_molec=`sed -n '8p' $3/$4/logs.txt`

echo "writing jobs"
python jobid_writer.py -pt $protein -fp $file_path -n_it $1 -jid $SLURM_JOB_NAME -jn $SLURM_JOB_NAME.txt

t_pos=$2    # total number of processers available
echo "Extracting labels"

if [ $sof = 'Glide' ]; then
   kw='r_i_docking_score'
elif [ $sof = 'FRED' ]; then
   kw='FRED Chemgauss4 score'
fi   

python scripts_2/extract_labels.py -n_it $1 -pt $protein -fp $file_path -t_pos $t_pos -score "$kw"

if [ $? != 0 ]; then
  echo "Extract_labels failed... terminating"
  exit
fi

part_gpu=$5

if [ $6 = $1 ]; then
   last='True'
else
   last='False'
fi

echo "Creating simple jobs"
python scripts_2/simple_job_models.py -n_it $1 -mdd $morgan_directory -time $time -file_path $file_path/$protein -nhp $nhp -titr $6 -n_mol $num_molec -pfm $7 -plm $8 -ct $rec -gp $part_gpu -tf_e $env -isl $last

cd $file_path/$protein/iteration_$1
rm model_no.txt
cd simple_job

echo "Running simple jobs"
#Executes all the files that were created in the simple_jobs directory
for f in *;do sbatch $f;done
```
The shell script logic is pretty straight forward, so I will move to the `extract_labels.py` script. 

#### `extract_labels.py`

```py
import argparse
import builtins as __builtin__
import glob
import gzip
import os
from contextlib import closing
from multiprocessing import Pool


# For debugging purposes only:
def print(*args, **kwargs):
    __builtin__.print('\t extract_L: ', end="")
    return __builtin__.print(*args, **kwargs)


parser = argparse.ArgumentParser()
parser.add_argument('-pt', '--project_name', required=True, help='Name of the DD project')
parser.add_argument('-fp', '--file_path', required=True, help='Path to the project directory, excluding project directory name')
parser.add_argument('-n_it', '--iteration_no', required=True, help='Number of current iteration')
parser.add_argument('-t_pos', '--tot_process', required=True, help='Number of CPUs to use for multiprocessing')
parser.add_argument('-score', '--score_keyword', required=True, help='Score keyword. Name of the field storing the docking score in the SDF files of docking results')

io_args = parser.parse_args()
n_it = int(io_args.iteration_no)
protein = io_args.project_name
file_path = io_args.file_path
tot_process = int(io_args.tot_process)
key_word = str(io_args.score_keyword)

# mol_key = 'ZINC'
print("Keyword: ", key_word)


def get_scores(ref):
    scores = []
    for line in ref:  # Looping through the molecules
        zinc_id = line.rstrip()
        line = ref.readline()
        # '$$$' signifies end of molecule info
        while line != '' and line[:4] != '$$$$':  # Looping through its information and saving scores

            tmp = line.rstrip().split('<')[-1]

            if key_word == tmp[:-1]:
                tmpp = float(ref.readline().rstrip())
                if tmpp > 50 or tmpp < -50:
                    print(zinc_id, tmpp)
                else:
                    scores.append([zinc_id, tmpp])

            line = ref.readline()
    return scores


def extract_glide_score(filen):
    scores = []
    try:
        # Opening the GNU compressed file
        with gzip.open(filen, 'rt') as ref:
            scores = get_scores(ref)

    except Exception as e:
        print('Exception occured in Extract_labels.py: ', e)
        # file is already decompressed
        with open(filen, 'r') as ref:
            scores = get_scores(ref)

    if 'test' in os.path.basename(filen):
        new_name = 'testing'
    elif 'valid' in os.path.basename(filen):
        new_name = 'validation'
    elif 'train' in os.path.basename(filen):
        new_name = 'training'
    else:
        print("Could not generate new training set")
        exit()

    with open(file_path + '/' + protein + '/iteration_' + str(n_it) + '/' + new_name + '_' + 'labels.txt', 'w') as ref:
        ref.write('r_i_docking_score' + ',' + 'ZINC_ID' + '\n')
        for z_id, gc in scores:
            ref.write(str(gc) + ',' + z_id + '\n')


if __name__ == '__main__':
    files = []
    iter_path = file_path + '/' + protein + '/iteration_' + str(n_it)

    # Checking to see if the labels have already been extracted:
    sets = ["training", "testing", "validation"]
    files_labels = glob.glob(iter_path + "/*_labels.txt")
    foundAll = True
    for s in sets:
        found = False
        print(s)
        for f in files_labels:
            set_name = f.split('/')[-1].split("_labels.txt")[0]
            if set_name == s:
                found = True
                print('Found')
                break
        if not found:
            foundAll = False
            print('Not Found')
            break
    if foundAll:
        print('Labels have already been extracted...')
        print('Remove "_labels.text" files in \"' + iter_path + '\" to re-extract')
        exit(0)

    path = iter_path + '/docked/*.sdf*'
    path_labels = iter_path + '/*labels*'

    for f in glob.glob(path):
        files.append(f)

    print("num files in", path, ":", len(files))
    print("Files:", [os.path.basename(f) for f in files])
    if len(files) == 0:
        print('NO FILES IN: ', path)
        print('CANCEL JOB...')
        exit(1)

    # Parallel running of the extract_glide_score() with each file path of the files array
    with closing(Pool(len(files))) as pool:
        pool.map(extract_glide_score, files)
```

##### Functions 

The `get_scores()` function is used to extract the scores from the SDF files. It is used *inside* `extract_glide_score()` like so: 

```py
scores = get_scores(ref)
```

It is defined like so: 
```py
def get_scores(ref):
    scores = []
    for line in ref:  # Looping through the molecules
        zinc_id = line.rstrip()
        line = ref.readline()
        # '$$$' signifies end of molecule info
        while line != '' and line[:4] != '$$$$':  # Looping through its information and saving scores
            tmp = line.rstrip().split('<')[-1]
            if key_word == tmp[:-1]:
                tmpp = float(ref.readline().rstrip())
                if tmpp > 50 or tmpp < -50:
                    print(zinc_id, tmpp)
                else:
                    scores.append([zinc_id, tmpp])
            line = ref.readline()
    return scores
```

This function will be called on an SDF file, we start with an empty scores list and start reading it like a text file. The first line is the molecule identifier (called `zinc_id` here). We strip the white space and save it. Then we read the next line. If the line is not empty and the first four characters are not `$$$$` (which signifies the end of the molecule information), we enter the while loop. We set tmp to the stripped line, split it on `<`, and take the last element. If the `key_word` (which is the score keyword) is equal to the `tmp` (minus the last character, which is `>`), *we read the next line* (`tmpp = float(ref.readline().rstrip())`), strip it, and convert it to a float. If the float is greater than 50 or less than -50, we print the `zinc_id` and the `tmpp` (which is the score). Otherwise, we append the `zinc_id` and the `tmpp` to the `scores` list. 

I think this is really clever, I always struggle with parsing SDF/text files for data like this. I think this is a really good way to do it. 

To better understand, I'll use an example sdf file for benzene. I am going to truncate the file slightly as we don't need all the data, we just need the structure of an sdf entry:

```sh
Benzene_Demo # Name and start of the molecule
                    3D
 Schrodinger Suite 2023-3.
 12 12  0  0  1  0            999 V2000
  -16.9063  -33.3752  -27.6873 C   0  0  0  0  0  0 # 3D coodrinates
  1  2  1  0  0  0 # Connectivity 
M  END
> <i_epik_Tot_Q> # Various property values
0

> <r_i_docking_score>
-4.60724452472058

$$$$ # End of molecule
```

When we start `for line in ref`, we are going to read this file line by line. The first line we see is going to be our molecule ID, and we save that. So we are here in the file: 

```
Benzene_Demo # Name and start of the molecule
```

Then, crucially, we advance the file pointer to the next line *which affects the for loop*. The for loop is now ahead one more increment. 

Now we are here: 

```
Benzene_Demo # Name and start of the molecule
                    3D
```

We check if the line is empty or if the first four characters are `$$$$`. They are not, so we enter the while loop.

We set `tmp` to the last element of the split line. This is going to be `r_i_docking_score>`. We check if `key_word` is equal to `tmp[:-1]`. If it is is we are going to read the line ahead of it to get the docking score (again we increment the for loop position) This is going to be `-4.60724452472058`. We are now here:

```
Benzene_Demo # Name and start of the molecule
                    3D
 Schrodinger Suite 2023-3.
 12 12  0  0  1  0            999 V2000
  -16.9063  -33.3752  -27.6873 C   0  0  0  0  0  0 # 3D coodrinates
  1  2  1  0  0  0 # Connectivity 
M  END
> <i_epik_Tot_Q> # Various property values
0

> <r_i_docking_score>
-4.60724452472058
```

 We strip it and convert it to a float. We check if it is greater than 50 or less than -50. It is not, so we append the `zinc_id` and the `tmpp` to the `scores` list. 

 Then we again increment the for loop pointer and we are at the `$$$$`:

 ```sh
Benzene_Demo # Name and start of the molecule
                    3D
 Schrodinger Suite 2023-3.
 12 12  0  0  1  0            999 V2000
  -16.9063  -33.3752  -27.6873 C   0  0  0  0  0  0 # 3D coodrinates
  1  2  1  0  0  0 # Connectivity 
M  END
> <i_epik_Tot_Q> # Various property values
0

> <r_i_docking_score>
-4.60724452472058

$$$$ # End of molecule
```

Now, we go back to the for loop, that will read the line *after* the `$$$$`. If we have more molecules, we will then get that ID. If we don't we have reached the end of the file.

When feeding the above explanation through coding LLMs, it doesn't like how I've phrased skipping lines/advancing the file pointer. However the practical outcome seems the same.

The rest of the script is actually pretty readable to me. We are essentially supplying paths so we can call `extract_scores()` on the desired files. We are also checking to see if the labels have already been extracted. If they have, we exit. I think I'd actually like this functionality in the other parts of the workflow as well, as a safe gaurd. 

We can now move on to the `simple_job_models.py` script. 

#### `simple_job_models.py`

The next portion of the shell script calls this line:

```sh
python scripts_2/simple_job_models.py -n_it $1 -mdd $morgan_directory -time $time -file_path $file_path/$protein -nhp $nhp -titr $6 -n_mol $num_molec -pfm $7 -plm $8 -ct $rec -gp $part_gpu -tf_e $env -isl $last
```

This is the `simple_job_models.py` script:

```py
import builtins as __builtin__
import pandas as pd
import numpy as np
import argparse
import glob
import time
import os

try:
    import __builtin__
except ImportError:
    # Python 3
    import builtins as __builtin__
    
# For debugging purposes only:
def print(*args, **kwargs):
    __builtin__.print('\t simple_jobs: ', end="")
    return __builtin__.print(*args, **kwargs)


START_TIME = time.time()


parser = argparse.ArgumentParser()
parser.add_argument('-n_it','--iteration_no',required=True,help='Number of current iteration')
parser.add_argument('-mdd','--morgan_directory',required=True,help='Path to the Morgan fingerprint directory for the database')
parser.add_argument('-time','--time',required=True,help='Time limit for training')
parser.add_argument('-file_path','--file_path',required=True,help='Path to the project directory, including project directory name')
parser.add_argument('-nhp','--number_of_hyp',required=True,help='Number of hyperparameters')
parser.add_argument('-titr','--total_iterations',required=True,help='Desired total number of iterations')

parser.add_argument('-isl','--is_last',required=True,help='True/False for is this last iteration')
parser.add_argument('-gp','--gpu_part',required=True,help='name(s) of GPU partitions')
parser.add_argument('-tf_e','--tensorflow_env',required=True,help='name of tensorflow-gpu environment')

# adding parameter for where to save all the data to:
parser.add_argument('-save', '--save_path', required=False, default=None)

# allowing for variable number of molecules to test and validate from:
parser.add_argument('-n_mol', '--number_mol', required=False, default=1000000, help='Size of validation/test set to be used')

parser.add_argument('-pfm', '--percent_first_mols', required=True, help='% of top scoring molecules to be considered as virtual hits in the first iteration (for standard DD run on 11 iterations, we recommend 1)')  # these two inputs must be percentages
parser.add_argument('-plm', '--percent_last_mols', required=True, help='% of top scoring molecules to be considered as virtual hits in the last iteration (for standard DD run on 11 iterations, we recommend 0.01)')


# Pass the threshold
parser.add_argument('-ct', '--recall', required=False, default=0.9, help='Recall, [0,1] range, default value 0.9')


funct_flags = parser.add_mutually_exclusive_group(required=False)
funct_flags.add_argument('-expdec', '--exponential_dec', required=False, default=-1) # must pass in the base number
funct_flags.add_argument('-polydec', '--polynomial_dec', required=False, default=-1) # You must also pass in to what power for this flag

io_args, extra_args = parser.parse_known_args()
n_it = int(io_args.iteration_no)
mdd = io_args.morgan_directory
time_model = io_args.time
nhp = int(io_args.number_of_hyp)
isl = io_args.is_last
titr = int(io_args.total_iterations)
rec = float(io_args.recall)
gpu_part = str(io_args.gpu_part)
env  = str(io_args.tensorflow_env)
protein = str(io_args.file_path).split('/')[-1]

num_molec = int(io_args.number_mol)

percent_first_mols = float(io_args.percent_first_mols)/100
percent_last_mols = float(io_args.percent_last_mols)/100

exponential_dec = int(io_args.exponential_dec)
polynomial_dec = int(io_args.polynomial_dec)

DATA_PATH = io_args.file_path   # Now == file_path/protein
SAVE_PATH = io_args.save_path
# if no save path is provided we just save it in the same location as the data
if SAVE_PATH is None: SAVE_PATH = DATA_PATH

# sums the first column and divides it by 1 million (this is our total database size)
t_mol = pd.read_csv(mdd+'/Mol_ct_file_%s.csv'%protein,header=None)[[0]].sum()[0]/1000000 # num of compounds in each file is mol_ct_file

cummulative = 0.25*n_it
dropout = [0.2, 0.5]
learn_rate = [0.0001]
bin_array = [2, 3]
wt = [2, 3]
if nhp < 144:
   bs = [256]
else:
    bs = [128, 256]
    
if nhp < 48:
    oss = [10]
elif nhp < 72:
    oss = [5, 10]
else:
    oss = [5, 10, 20]

if nhp < 24:
    num_units = [1500, 2000]
else:
    num_units = [100, 1500, 2000]

try:
    os.mkdir(SAVE_PATH+'/iteration_'+str(n_it)+'/simple_job')
except OSError: # catching file exists error
    pass

# Clearing up space from previous iteration
for f in glob.glob(SAVE_PATH+'/iteration_'+str(n_it)+'/simple_job/*'):
    os.remove(f)

scores_val = []
with open(DATA_PATH+'/iteration_'+str(1)+'/validation_labels.txt','r') as ref:
    ref.readline()  # first line is ignored
    for line in ref:
        scores_val.append(float(line.rstrip().split(',')[0]))

scores_val = np.array(scores_val)

first_mols = int(100*t_mol/13) if percent_first_mols == -1.0 else int(percent_first_mols * len(scores_val))

print(first_mols)

last_mols = 100 if percent_last_mols == -1.0 else int(percent_last_mols * len(scores_val))

print(last_mols)

if n_it==1:
    # 'good_mol' is the number of top scoring molecules to save at the end of the iteration
    good_mol = first_mols
else:
    if exponential_dec != -1:
        good_mol = int() #TODO: create functions for these
    elif polynomial_dec != -1:
        good_mol = int()
    else:
        good_mol = int(((last_mols-first_mols)*n_it + titr*first_mols-last_mols)/(titr-1))     # linear decrease as interations increase


print(isl)
#If this is the last iteration then we save only 100 molecules
if isl == 'True':
    good_mol = 100 if percent_last_mols == -1.0 else int(percent_last_mols * len(scores_val))

cf_start = np.mean(scores_val)  # the mean of all the docking scores (labels) of the validation set:
t_good = len(scores_val)

# we decrease the threshold value until we have our desired num of mol left.
while t_good > good_mol: 
    cf_start -= 0.005
    t_good = len(scores_val[scores_val<cf_start])

print('Threshold (cutoff):',cf_start)
print('Molec under threshold:', t_good)
print('Goal molec:', good_mol)
print('Total molec:', len(scores_val))

all_hyperparas = []

for o in oss:   # Over Sample Size
    for batch in bs:
        for nu in num_units:
            for do in dropout:
                for lr in learn_rate:
                    for ba in bin_array:
                        for w in wt:    # Weight
                            all_hyperparas.append([o,batch,nu,do,lr,ba,w,cf_start])

print('Total hyp:', len(all_hyperparas))

# Creating all the jobs for each hyperparameter combination:

other_args = ' '.join(extra_args) + '-rec {} -n_it {} -t_mol {} --data_path {} --save_path {} -n_mol {}'.format(rec, n_it, t_mol, DATA_PATH, SAVE_PATH, num_molec)
print(other_args)
count = 1
for i in range(len(all_hyperparas)):
    with open(SAVE_PATH+'/iteration_'+str(n_it)+'/simple_job/simple_job_'+str(count)+'.sh', 'w') as ref:
        ref.write('#!/bin/bash\n')
        ref.write('#SBATCH --ntasks=1\n')
        ref.write('#SBATCH --gres=gpu:1\n')
        ref.write('#SBATCH --cpus-per-task=1\n')
        ref.write('#SBATCH --job-name=phase_4\n')
        ref.write('#SBATCH --mem=0               # memory per node\n')
        ref.write('#SBATCH --partition=%s\n'%gpu_part)
        ref.write('#SBATCH --time='+time_model+'            # time (DD-HH:MM)\n')
        ref.write('\n')
        cwd = os.getcwd()
        ref.write('cd {}/scripts_2\n'.format(cwd))
        hyp_args = '-os {} -bs {} -num_units {} -dropout {} -learn_rate {} -bin_array {} -wt {} -cf {}'.format(*all_hyperparas[i])
        ref.write('source ~/.bashrc\n')
        ref.write('conda activate %s\n'%env)
        ref.write('python -u progressive_docking.py ' + hyp_args + ' ' + other_args)
        ref.write("\n echo complete")
    count += 1
    
print('Runtime:', time.time() - START_TIME)
```

A fairly long script, but not incredibly difficult. The goal of this script is to write out the `progressive_docking.py` script with all the hyperparameters we want to test. I think somewhere here is where we would integrate something like `Optuna` if we wanted to. Let's go through the script.

I will skip over the imports and the `print` function. 

The arguments are mostly what we expect, but there are a few new ones - particularly `-isl` or `--is_last` which is a boolean flag that changes behavior on whether we are at the last iteration or not. 

Additionally, we start seeing terms that are related to the machine learning aspect of this workflow. Two terms I am not familiar with are:

```py
funct_flags = parser.add_mutually_exclusive_group(required=False)
funct_flags.add_argument('-expdec', '--exponential_dec', required=False, default=-1) # must pass in the base number
funct_flags.add_argument('-polydec', '--polynomial_dec', required=False, default=-1) # You must also pass in to what power for this flag
```

So first, we have an argparse functionality I wasn't aware of, the ability to add mutually exclusive arguments. Only one of these can be set. 

However, this doesn't even seem to matter because the script never sets these terms when it is called earlier in the shell script: 

```sh
python scripts_2/simple_job_models.py -n_it $1 -mdd $morgan_directory -time $time -file_path $file_path/$protein -nhp $nhp -titr $6 -n_mol $num_molec -pfm $7 -plm $8 -ct $rec -gp $part_gpu -tf_e $env -isl $last
```
Here are the terms we do call and their definitions:

- `-n_it` `$1`: Number of the current iteration (passed from the first positional argument in the shell script).
- `-mdd` `$morgan_directory`: Path to the directory containing Morgan fingerprints.
- `-time` `$time`: Time limit for the training process.
- `-file_path` `$file_path/$protein`: Path to the project directory, including the project directory name, with a specific protein.
- `-nhp` `$nhp`: Number of hyperparameters.
- `-titr` `$6`: Total desired number of iterations (passed from the sixth positional argument).
- `-n_mol` `$num_molec`: Number of molecules to consider for validation/testing.
- `-pfm` `$7`: Percentage of top-scoring molecules to be considered as virtual hits in the first iteration (passed from the seventh positional argument).
- `-plm` `$8`: Percentage of top-scoring molecules to be considered as virtual hits in the last iteration (passed from the eighth positional argument).
- `-ct` `$rec`: Recall, a threshold value in the range [0,1].
- `-gp` `$part_gpu`: Name(s) of GPU partitions.
- `-tf_e` `$env`: Name of the TensorFlow-GPU environment.
- `-isl` `$last`: A flag indicating whether this is the last iteration.

Reading the code more closely regarding the mutually exclusive args, I see this part, which seems to indicate these functions are not yet implemented: 

```py 
if n_it==1:
    good_mol = first_mols
else:
    if exponential_dec != -1:
        good_mol = int() #TODO: create functions for these
    elif polynomial_dec != -1:
        good_mol = int()
    else:
        good_mol = int(((last_mols-first_mols)*n_it + titr*first_mols-last_mols)/(titr-1)) # linear decrease as iterations increase
```

So we will hit the `else:` condition. Here, `good_mol` is set to a value that is a linear function of the number of iterations. 

I think this governs the fact that the amount of data we are considering decreases as we go through more iterations, because we need to challenge the model with more difficult data to make it better. 

This something we probably need to consider during the strain incorporation process. 

Here are some of the args that are setting `good_mol`:

```py
first_mols = int(100*t_mol/13) if percent_first_mols == -1.0 else int(percent_first_mols * len(scores_val))

print(first_mols)

last_mols = 100 if percent_last_mols == -1.0 else int(percent_last_mols * len(scores_val))
```

I have never set `percent_first_mols` to -1, so I am not sure exactly why that triggers, but generally first_mols is going to be in the `else` statement. 

`percent_first_mols` and `percent_last_mols` are set in the `argparse` portion of the script. They are defined to be percentages, i.e. the percentage of top scoring molecules with in the bounds of your first and last iterations. You may recall from the original paper that with 10 iterations we move from 1% to 0.01% percent from iteration 1 -> iteration 11. These are in the unit type of *percentages*. Here's the `argparse` portion: 

```py
parser.add_argument('-pfm', '--percent_first_mols', required=True, help='% of top scoring molecules to be considered as virtual hits in the first iteration (for standard DD run on 11 iterations, we recommend 1)')  # these two inputs must be percentages
parser.add_argument('-plm', '--percent_last_mols', required=True, help='% of top scoring molecules to be considered as virtual hits in the last iteration (for standard DD run on 11 iterations, we recommend 0.01)')
```

So, with that in mind - let's think about the "default" parameters here. We can simplify that statement: 

`first_mols = int(1 * len(scores_val))`

`last_mols = int(0.01 * len(scores_val))`

Now we need `scores_val`. 

That is set a little earlier in the script, it uses the validation_labels.txt in iteration 1 and does not change. 

```py
scores_val = []
with open(DATA_PATH+'/iteration_'+str(1)+'/validation_labels.txt','r') as ref:
    ref.readline()  # first line is ignored
    for line in ref:
        scores_val.append(float(line.rstrip().split(',')[0]))

scores_val = np.array(scores_val)
```
Let's again assume default parameters and say we have 1,000,000 (1M) molecules in our validation set. The first line is assumed to be the header, so we can ignore it. 

Now, we can simplify again and do the arthimetic: 

`first_mols = int(1 * 1,000,000) = 1,000,000`
`last_mols = int(0.01 * 1,000,000) = 10,000`

Now, let's return to `good_mol`:

```py
good_mol = int(((last_mols-first_mols)*n_it + titr*first_mols-last_mols)/(titr-1)) # linear decrease as iterations increase
```

Let's assume we are on iteration 1 and substitute in our values:

`good_mol = int(((10,000 - 1,000,000) * 1 + 11 * 1,000,000 - 10,000)/(11-1))`

This simplifies down to (I think, just ran it through WolframAlpha): 1,000,000 (1M)

However, this never actually runs when `n_it == 1` because of the `if` statement we saw earlier!

So, let's assume we are on iteration 2. 

`good_mol = int(((10,000 - 1,000,000) * 2 + 11 * 1,000,000 - 10,000)/(11-1))`

This simplifies down to: 901,000. This is 90.1% of the original data set, but I did expect it to be an even 90%. Perhaps this will make sense later. Let's quickly make a table of n_it and good_mol:

| n_it | good_mol |
|------|----------|
| 1    | 1,000,000 |
| 2    | 901,000 |
| 3    | 802,000 |
| 4    | 703,000 |
| 5    | 604,000 |
| 6    | 505,000 |
| 7    | 406,000 |
| 8    | 307,000 |
| 9    | 208,000 |
| 10?   | 109,000? |

We also have a special condition regarding if it is the last iteration. The syntax/word choice here I am using is probably a little wrong, as I I feel like by n_it == 10, we are on the last iteration. So it will trigger there. 

So I think it will be: 

| n_it | good_mol |
|------|----------|
| 10   | 10,000   |

```py
print(isl)
#If this is the last iteration then we save only 100 molecules
if isl == 'True':
    good_mol = 100 if percent_last_mols == -1.0 else int(percent_last_mols * len(scores_val))
```

The comment is probably wrong, as you can set `percent_last_mols` to a number other than 0.01 and your `scores_val` maybe a different value than 1,000,000. But either way, the logic does make sense, so let's continue on. 

Something I just noticed and forgot however is this: 

```py
percent_first_mols = float(io_args.percent_first_mols)/100
percent_last_mols = float(io_args.percent_last_mols)/100
```
So, we are dividing by 100. So, let's update our table (keep in mind n_it == 1 is a special case):

| n_it | good_mol |
|------|----------|
| 1    | 1,000,000 |
| 2    | 9,010 |
| 3    | 8,020 |
| 4    | 7,030 |
| 5    | 6,040 |
| 6    | 5,050 |
| 7    | 4,060 |
| 8    | 3,070 |
| 9    | 2,080 |
| 10?   | 1,090? |
| 10 (isl=True)  | 100 |

I'll have to rewrite this section later, I just don't want to get hung up on that now. But because `percent_first/last_mols` are divided by 100, the earlier arithmetic that governs `first/last_mol` should have been divided by 100:

```py
first_mols = int(100*t_mol/13) if percent_first_mols == -1.0 else int(percent_first_mols * len(scores_val))

print(first_mols)

last_mols = 100 if percent_last_mols == -1.0 else int(percent_last_mols * len(scores_val))
```

So here it would have been: 

`first_mols = int(0.01 * len(scores_val)) = 10,000`
`last_mols = int(0.0001 * len(scores_val)) = 100`

Then, when we get to `good_mol` we would have: 

`good_mol = int(((100 - 10,000) * 1 + 11 * 10,000 - 100)/(11-1)) = 10,000`

Again, this doesn't actually trigger when `n_it == 1` because of the `if` statement. So let's evaluate when `n_it == 2`:

`good_mol = int(((100 - 10,000) * 2 + 11 * 10,000 - 100)/(11-1)) = 9,010`

Now we get to the next block of code: 

```py
cf_start = np.mean(scores_val)  # the mean of all the docking scores (labels) of the validation set:
t_good = len(scores_val)
```

Pretty straightforward.

Now the next chunk of code is a `while` loop. 
```py
# we decrease the threshold value until we have our desired num of mol left.
while t_good > good_mol: 
    cf_start -= 0.005
    t_good = len(scores_val[scores_val<cf_start])
```
New operator I haven't seen: 

`cf_start -= 0.005` is equivalent to: 
`cf_start = cf_start - 0.005`

So remember again that `scores_val` is a constant. So, let's walk through this loop. 

In iteration 1, we have: 

`t_good` = 1,000,000
`good_mol` = 1,000,000

Right now, `t_good` is *not* greater than `good_mol` so we do not enter the loop. 

In iteration 2, we have:

`t_good` = 1,000,000
`good_mol` = 9,010

Now, `t_good` is greater than `good_mol` so we enter the loop. 

`cf_start` is set to the mean of `scores_val`, let's just pretend that the average is `-5.0` for now. 

So we will subtract `0.005` from `-5.0` and get `-5.005` (a better docking score). 

Then we will set `t_good` to the number of elements in `scores_val` that are less than `-5.005`. 

So, let's say that there are 9,000 elements in `scores_val` that are less than `-5.005`. 

So, we will set `t_good` to 9,000 as exit the while loop and print the results:

```py
print('Threshold (cutoff):',cf_start)
print('Molec under threshold:', t_good)
print('Goal molec:', good_mol)
print('Total molec:', len(scores_val))
```

So, we will get: 

```
Threshold (cutoff): -5.005
Molec under threshold: 9,000
Goal molec: 9,010
Total molec: 1,000,000
```
I am a little confused on t_good versus good_mol, but maybe the word "goal" is tripping me up. 

This again stresses to me that these scripts should be outputting logs associated with them and saved to an appropriate place, where the logs are these print statements (and the command that was executed). 

Now we get to the hyperparameter section. 

```py
all_hyperparas = []

for o in oss:   # Over Sample Size
    for batch in bs:
        for nu in num_units:
            for do in dropout:
                for lr in learn_rate:
                    for ba in bin_array:
                        for w in wt:    # Weight
                            all_hyperparas.append([o,batch,nu,do,lr,ba,w,cf_start])

print('Total hyp:', len(all_hyperparas))
```
I find this for loop pretty crazy but I think I understand at least what its trying to do. We are going to be appending a list of lists to `all_hyperparas` with all the combinations we can make within the for loops. 

Next we have a pretty useful code snippet I'd like to remember. We are about to write out the hyperparameter combinations to a series of files that will run `progressive_docking.py` with them. But, `progressive_docking.py` also needs a set of arguments that don't change, so they are defined here: 

```py
other_args = ' '.join(extra_args) + '-rec {} -n_it {} -t_mol {} --data_path {} --save_path {} -n_mol {}'.format(rec, n_it, t_mol, DATA_PATH, SAVE_PATH, num_molec)
print(other_args)
```
This will print out a string that looks like this:

`-pfm 0.01 -plm 0.0001 -rec 0.5 -n_it 1 -t_mol 1000000 --data_path data/path/ --save_path also/data/path -n_mol 1000000`

Then, we will use that when we get to the final section: 

```py
count = 1
for i in range(len(all_hyperparas)):
    with open(SAVE_PATH+'/iteration_'+str(n_it)+'/simple_job/simple_job_'+str(count)+'.sh', 'w') as ref:
        ref.write('#!/bin/bash\n')
        ref.write('#SBATCH --ntasks=1\n')
        ref.write('#SBATCH --gres=gpu:1\n')
        ref.write('#SBATCH --cpus-per-task=1\n')
        ref.write('#SBATCH --job-name=phase_4\n')
        ref.write('#SBATCH --mem=0               # memory per node\n')
        ref.write('#SBATCH --partition=%s\n'%gpu_part)
        ref.write('#SBATCH --time='+time_model+'            # time (DD-HH:MM)\n')
        ref.write('\n')
        cwd = os.getcwd()
        ref.write('cd {}/scripts_2\n'.format(cwd))
        hyp_args = '-os {} -bs {} -num_units {} -dropout {} -learn_rate {} -bin_array {} -wt {} -cf {}'.format(*all_hyperparas[i])
        ref.write('source ~/.bashrc\n')
        ref.write('conda activate %s\n'%env)
        ref.write('python -u progressive_docking.py ' + hyp_args + ' ' + other_args)
        ref.write("\n echo complete")
    count += 1
    
print('Runtime:', time.time() - START_TIME)
```

Which get's us to the portion where we can now do the machine learning related stuff!