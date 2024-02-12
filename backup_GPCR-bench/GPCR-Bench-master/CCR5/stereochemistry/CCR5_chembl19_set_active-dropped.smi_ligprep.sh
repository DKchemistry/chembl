#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME CCR5_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi CCR5_chembl19_set_active-dropped.smi -osd CCR5_chembl19_set_active-dropped.smi.sdf
