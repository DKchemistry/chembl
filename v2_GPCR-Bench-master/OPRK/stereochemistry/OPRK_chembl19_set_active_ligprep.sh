#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME OPRK_chembl19_set_active_ligprep -HOST localhost:1 -ismi OPRK_chembl19_set_active_sc.smi -osd OPRK_chembl19_set_active.sdf
