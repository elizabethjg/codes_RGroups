#!/bin/bash
#SBATCH --job-name=Mbin1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /home/elizabeth/codes_RGroups

export OMP_NUM_THREADS=10

module load softwares/anaconda2/5.0-intel-2018.0
source activate py2env

srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
-file 'profile_Mbin1.fits' -ncores 15 -Yanmiss 'True'
srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
-file 'profile_Mbin1_conM.fits' -ncores 15 -Yanmiss 'True'
