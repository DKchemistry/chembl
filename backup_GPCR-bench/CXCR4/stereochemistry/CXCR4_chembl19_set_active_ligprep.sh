#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME CXCR4_chembl19_set_active_ligprep -HOST localhost:1 -ismi CXCR4_chembl19_set_active_sc.smi -osd CXCR4_chembl19_set_active.sdf
