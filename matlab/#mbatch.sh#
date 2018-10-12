#!/bin/bash

outfile=/project/ecog/emily/PF/$JOB_ID
if ! [ -x $outfile ]; then mkdir $outfile; fi
matlab -nodisplay -singleCompThread -r "outfile = '$outfile/$SGE_TASK_ID', particle_filter, exit"

