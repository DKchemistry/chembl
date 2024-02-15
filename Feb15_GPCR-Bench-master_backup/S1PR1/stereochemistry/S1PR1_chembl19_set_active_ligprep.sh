#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME S1PR1_chembl19_set_active_ligprep -HOST localhost:1 -ismi S1PR1_chembl19_set_active_sc.smi -osd S1PR1_chembl19_set_active.sdf
