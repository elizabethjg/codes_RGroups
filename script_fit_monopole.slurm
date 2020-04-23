#!/bin/bash
#SBATCH --job-name=N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

source activate py2env

srun -u python fit_profile_monopole_misscentred.py -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' -file 'profile_N1.fits' -ncores 10
 