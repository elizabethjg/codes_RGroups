#!/bin/bash
#SBATCH --job-name=Mbin1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=2-00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /home/elizabeth/codes_RGroups

export OMP_NUM_THREADS=10

module load softwares/anaconda2/5.0-intel-2018.0
source activate py2env

srun -u python fit_profile_monopole_misscentred_pcc2h.py -Mbin $SLURM_JOB_NAME

##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_indcat/' \
##-file 'profile_'$SLURM_JOB_NAME'.fits' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_indcat/' \
##-file 'profile_'$SLURM_JOB_NAME'.fits' -Yanmiss 'True' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_new/' \
##-file 'profile_'$SLURM_JOB_NAME'_cM.fits' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_newanalysis/' \
##-file 'profile_'$SLURM_JOB_NAME'.fits' -Yanmiss 'True' -ncores 15                                    
##srun -u python fit_profile_monopole_misscentred_pcc.py -folder '/home/elizabeth/profiles_new/' \
##-file 'profile_'$SLURM_JOB_NAME'_cM.fits' -Yanmiss 'True' -ncores 15                                    
