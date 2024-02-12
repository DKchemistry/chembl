#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME OPRM_chembl19_set_active_ligprep -HOST localhost:1 -ismi OPRM_chembl19_set_active_sc.smi -osd OPRM_chembl19_set_active.sdf
