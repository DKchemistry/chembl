#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME PAR1_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi PAR1_chembl19_set_decoy-dropped.smi -osd PAR1_chembl19_set_decoy-dropped.smi.sdf
