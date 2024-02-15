#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME OX2R_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi OX2R_chembl19_set_active-dropped.smi -osd OX2R_chembl19_set_active-dropped.smi.sdf
