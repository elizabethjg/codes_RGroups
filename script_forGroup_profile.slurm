#!/bin/bash
#SBATCH --job-name=profileN1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=56
export MKL_NUM_THREADS=56

source activate py2env

srun python -u forGroup_profile.py -sample 'N1' -N_min 1 -N_max 2 -ncores 56
 
