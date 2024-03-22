#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME CXCR4_chembl19_set_decoy_ligprep -HOST localhost:5 -ismi CXCR4_chembl19_set_decoy_sc.smi -osd CXCR4_chembl19_set_decoy.sdf
