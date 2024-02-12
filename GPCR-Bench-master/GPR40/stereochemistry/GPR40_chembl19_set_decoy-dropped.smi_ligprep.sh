#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME GPR40_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi GPR40_chembl19_set_decoy-dropped.smi -osd GPR40_chembl19_set_decoy-dropped.smi.sdf
