#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME ACM2_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi ACM2_chembl19_set_active-dropped.smi -osd ACM2_chembl19_set_active-dropped.smi.sdf
