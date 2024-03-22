#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME MGLUR1_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi MGLUR1_chembl19_set_decoy_sc.smi -osd MGLUR1_chembl19_set_decoy.sdf
