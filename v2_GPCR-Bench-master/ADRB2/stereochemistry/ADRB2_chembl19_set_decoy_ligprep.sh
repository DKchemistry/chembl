#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME ADRB2_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi ADRB2_chembl19_set_decoy_sc.smi -osd ADRB2_chembl19_set_decoy.sdf
