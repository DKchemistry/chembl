#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME ADRB1_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi ADRB1_chembl19_set_active-dropped.smi -osd ADRB1_chembl19_set_active-dropped.smi.sdf
