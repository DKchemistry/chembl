#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME HRH1_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi HRH1_chembl19_set_decoy_sc.smi -osd HRH1_chembl19_set_decoy.sdf
