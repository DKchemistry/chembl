#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME OPRK_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi OPRK_chembl19_set_decoy-dropped.smi -osd OPRK_chembl19_set_decoy-dropped.smi.sdf
