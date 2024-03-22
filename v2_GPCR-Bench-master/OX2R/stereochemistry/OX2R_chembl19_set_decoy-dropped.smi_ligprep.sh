#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME OX2R_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi OX2R_chembl19_set_decoy-dropped.smi -osd OX2R_chembl19_set_decoy-dropped.smi.sdf
