#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME CCR5_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi CCR5_chembl19_set_decoy-dropped.smi -osd CCR5_chembl19_set_decoy-dropped.smi.sdf
