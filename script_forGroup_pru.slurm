#!/bin/bash
#SBATCH --job-name=Mbin_ncm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=6:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=30
export MKL_NUM_THREADS=30

source activate py2env
srun python -u forGroup_profile.py -sample 'N1_Mbin6'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 30
