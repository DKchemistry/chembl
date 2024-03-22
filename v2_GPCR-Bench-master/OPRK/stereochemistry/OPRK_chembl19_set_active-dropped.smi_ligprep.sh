#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME OPRK_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi OPRK_chembl19_set_active-dropped.smi -osd OPRK_chembl19_set_active-dropped.smi.sdf
