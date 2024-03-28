#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME S1PR1_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi S1PR1_chembl19_set_decoy-dropped.smi -osd S1PR1_chembl19_set_decoy-dropped.smi.sdf
