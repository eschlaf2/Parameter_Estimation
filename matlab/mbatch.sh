#!/bin/bash

matlab -nodisplay -singleCompThread -r "outfile = '/project/ecog/emily/PF/$SGE_TASK_ID', particle_filter, save(outfile, sim, estimates), exit"