#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME AA2AR_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi AA2AR_chembl19_set_decoy_sc.smi -osd AA2AR_chembl19_set_decoy.sdf
