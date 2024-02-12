#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME DRD3_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi DRD3_chembl19_set_active-dropped.smi -osd DRD3_chembl19_set_active-dropped.smi.sdf
