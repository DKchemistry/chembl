#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME MGLUR5_chembl19_set_active_ligprep -HOST localhost:1 -ismi MGLUR5_chembl19_set_active_sc.smi -osd MGLUR5_chembl19_set_active.sdf
