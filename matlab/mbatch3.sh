#!/bin/bash

outfile=/projectnb/ecog/emily/PF/$JOB_ID
if ! [ -x $outfile ]; then mkdir $outfile; fi
matlab -nodisplay -r "outfile = '$outfile/$SGE_TASK_ID', pf_settings3, particle_filter, exit"

