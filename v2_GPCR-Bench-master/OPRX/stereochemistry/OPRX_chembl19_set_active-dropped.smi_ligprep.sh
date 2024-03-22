#!/bin/bash

$SCHRODINGER/ligprep -inp ligprep_noepik.inp -NJOBS 1 -JOBNAME OPRX_chembl19_set_active-dropped.smi_ligprep -HOST localhost:1 -ismi OPRX_chembl19_set_active-dropped.smi -osd OPRX_chembl19_set_active-dropped.smi.sdf
