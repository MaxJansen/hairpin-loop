#!/bin/bash
# Finds the path and filename for all hairpins in a directory:
# Run the code in a 'site_0' directory, for example (for 5GGS_Z):
# /work/lpdi/users/anthony/masif_runs/masif/masif_seed_search/data/masif_targets/targets/5GGS_Z/out_hairpins/5GGS_Z/site_0
$TARGET=$(echo $1 | rev | cut -d'/' -f 2 | rev)
find $1 -type f -name "*.pdb" > ${TARGET}_selected_seeds.txt
