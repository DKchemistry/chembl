#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME DRD3_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi DRD3_chembl19_set_decoy_sc.smi -osd DRD3_chembl19_set_decoy.sdf
