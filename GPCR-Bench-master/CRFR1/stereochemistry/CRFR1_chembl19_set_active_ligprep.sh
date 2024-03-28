#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME CRFR1_chembl19_set_active_ligprep -HOST localhost:1 -ismi CRFR1_chembl19_set_active_sc.smi -osd CRFR1_chembl19_set_active.sdf
