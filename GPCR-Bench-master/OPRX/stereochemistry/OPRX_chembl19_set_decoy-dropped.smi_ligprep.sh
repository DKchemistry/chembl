#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME OPRX_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi OPRX_chembl19_set_decoy-dropped.smi -osd OPRX_chembl19_set_decoy-dropped.smi.sdf
