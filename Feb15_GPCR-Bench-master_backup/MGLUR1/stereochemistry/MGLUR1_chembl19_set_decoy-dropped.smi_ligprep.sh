#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME MGLUR1_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi MGLUR1_chembl19_set_decoy-dropped.smi -osd MGLUR1_chembl19_set_decoy-dropped.smi.sdf
