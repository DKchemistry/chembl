#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME 5HT2B_chembl19_set_active_ligprep -HOST localhost:1 -ismi 5HT2B_chembl19_set_active_sc.smi -osd 5HT2B_chembl19_set_active.sdf
