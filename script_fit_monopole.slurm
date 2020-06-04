#!/bin/bash
#SBATCH --job-name=f1_Mbin6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --time=16:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

source activate py2env

srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
-file 'profile_Mbin1.fits' -ncores 15 -Yanmiss 'True'
srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
-file 'profile_Mbin1_conM.fits' -ncores 15 -Yanmiss 'True'
