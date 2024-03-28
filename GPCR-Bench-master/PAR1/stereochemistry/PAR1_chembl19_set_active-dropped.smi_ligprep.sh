#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME PAR1_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi PAR1_chembl19_set_active-dropped.smi -osd PAR1_chembl19_set_active-dropped.smi.sdf
