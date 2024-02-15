#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 5 -JOBNAME ACM3_chembl19_set_decoy-dropped.smi_ligprep -HOST localhost:5 -ismi ACM3_chembl19_set_decoy-dropped.smi -osd ACM3_chembl19_set_decoy-dropped.smi.sdf
