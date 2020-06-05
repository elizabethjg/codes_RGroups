#!/bin/bash
cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/
source activate py2env
python -u forGroup_profile.py -sample 'N1_Mbin6'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 30
