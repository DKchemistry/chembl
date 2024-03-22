#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME P2Y12_chembl19_set_active_ligprep -HOST localhost:1 -ismi P2Y12_chembl19_set_active_sc.smi -osd P2Y12_chembl19_set_active.sdf
