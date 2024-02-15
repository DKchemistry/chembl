#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME DRD3_chembl19_set_active_ligprep -HOST localhost:1 -ismi DRD3_chembl19_set_active_sc.smi -osd DRD3_chembl19_set_active.sdf
