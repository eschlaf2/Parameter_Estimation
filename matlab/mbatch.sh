#!/bin/bash

matlab -nodisplay -singleCompThread -r "outfile = '/project/ecog/emily/PF/$1_$SGE_TASK_ID', particle_filter, save(outfile, sim, estimates), exit"
