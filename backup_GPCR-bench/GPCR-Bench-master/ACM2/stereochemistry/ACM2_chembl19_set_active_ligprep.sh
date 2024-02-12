#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME ACM2_chembl19_set_active_ligprep -HOST localhost:1 -ismi ACM2_chembl19_set_active_sc.smi -osd ACM2_chembl19_set_active.sdf
