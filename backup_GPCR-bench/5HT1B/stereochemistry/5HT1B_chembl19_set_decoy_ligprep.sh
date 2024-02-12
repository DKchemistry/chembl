#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME 5HT1B_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi 5HT1B_chembl19_set_decoy_sc.smi -osd 5HT1B_chembl19_set_decoy.sdf
