#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME OPRK_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi OPRK_chembl19_set_decoy_sc.smi -osd OPRK_chembl19_set_decoy.sdf
