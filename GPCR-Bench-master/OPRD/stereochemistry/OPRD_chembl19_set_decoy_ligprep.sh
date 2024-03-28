#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME OPRD_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi OPRD_chembl19_set_decoy_sc.smi -osd OPRD_chembl19_set_decoy.sdf
