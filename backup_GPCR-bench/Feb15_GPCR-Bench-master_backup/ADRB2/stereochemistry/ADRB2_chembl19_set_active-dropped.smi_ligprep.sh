#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME ADRB2_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi ADRB2_chembl19_set_active-dropped.smi -osd ADRB2_chembl19_set_active-dropped.smi.sdf
