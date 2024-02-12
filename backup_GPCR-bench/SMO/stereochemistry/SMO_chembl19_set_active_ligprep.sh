#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME SMO_chembl19_set_active_ligprep -HOST localhost:1 -ismi SMO_chembl19_set_active_sc.smi -osd SMO_chembl19_set_active.sdf
