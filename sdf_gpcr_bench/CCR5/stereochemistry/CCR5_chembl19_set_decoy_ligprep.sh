#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME CCR5_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi CCR5_chembl19_set_decoy_sc.smi -osd CCR5_chembl19_set_decoy.sdf
