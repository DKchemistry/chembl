#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME ADRB1_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi ADRB1_chembl19_set_decoy-dropped.smi -osd ADRB1_chembl19_set_decoy-dropped.smi.sdf
