#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME HRH1_chembl19_set_active_ligprep -HOST localhost:1 -ismi HRH1_chembl19_set_active_sc.smi -osd HRH1_chembl19_set_active.sdf
