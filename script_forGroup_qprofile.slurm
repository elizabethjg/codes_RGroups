#!/bin/bash
#SBATCH --job-name=qN1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=46
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=46
export MKL_NUM_THREADS=46

source activate py2env

srun python -u forGroup_qprofile.py -sample 'N1qBG' -N_min 1 -N_max 2 -C_BG_min 2.73 -ROUT 1000. -ncores 46
 
