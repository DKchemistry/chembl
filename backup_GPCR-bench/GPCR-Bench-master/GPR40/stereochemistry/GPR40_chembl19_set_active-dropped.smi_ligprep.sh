#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME GPR40_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi GPR40_chembl19_set_active-dropped.smi -osd GPR40_chembl19_set_active-dropped.smi.sdf
