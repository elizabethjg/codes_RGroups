#!/bin/bash
#SBATCH --job-name=qMbin
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

srun python -u forGroup_qprofile.py -sample 'Mbin125-14' -lMH_min 12.5 -lMH_max 14.0  -C_BG_min 0. -C_BG_max 100. -ncores 56
 
