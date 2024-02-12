#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME P2Y12_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi P2Y12_chembl19_set_decoy-dropped.smi -osd P2Y12_chembl19_set_decoy-dropped.smi.sdf
