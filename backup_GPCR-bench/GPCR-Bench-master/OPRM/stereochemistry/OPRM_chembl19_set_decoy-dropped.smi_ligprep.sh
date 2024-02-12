#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME OPRM_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi OPRM_chembl19_set_decoy-dropped.smi -osd OPRM_chembl19_set_decoy-dropped.smi.sdf
