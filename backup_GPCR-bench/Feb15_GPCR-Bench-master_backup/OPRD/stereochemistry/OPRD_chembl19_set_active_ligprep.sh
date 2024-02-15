#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME OPRD_chembl19_set_active_ligprep -HOST localhost:1 -ismi OPRD_chembl19_set_active_sc.smi -osd OPRD_chembl19_set_active.sdf
