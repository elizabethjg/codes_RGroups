#!/bin/bash
#SBATCH --job-name=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

. /etc/profile
source $HOME/.bashrc

cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

conda activate py2env

srun -u python fit_profile_monopole_misscentred_pcc2h.py -Mbin $SLURM_JOB_NAME

##srun -u python fit_profile_monopole_misscentred_pcc.py -file 'profile_'$SLURM_JOB_NAME'.fits' 

##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
-file 'profile_'$SLURM_JOB_NAME'.fits' -Yanmiss 'True' -ncores 15                                    
