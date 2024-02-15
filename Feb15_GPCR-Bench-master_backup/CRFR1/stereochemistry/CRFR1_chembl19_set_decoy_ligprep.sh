#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME CRFR1_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi CRFR1_chembl19_set_decoy_sc.smi -osd CRFR1_chembl19_set_decoy.sdf
