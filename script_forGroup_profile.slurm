#!/bin/bash
#SBATCH --job-name=CS82
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=49
#SBATCH --time=12:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=49
export MKL_NUM_THREADS=49

conda activate py2env


srun python -u forGroup_profile_CS82.py -sample 'Mbin1'$SLURM_JOB_NAME  -lMH_min 12.5 -lMH_max 12.7 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin2'$SLURM_JOB_NAME  -lMH_min 12.7 -lMH_max 12.9 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin3'$SLURM_JOB_NAME  -lMH_min 12.9 -lMH_max 13.1 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin4'$SLURM_JOB_NAME  -lMH_min 13.1 -lMH_max 13.3 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin5'$SLURM_JOB_NAME  -lMH_min 13.3 -lMH_max 13.5 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin6'$SLURM_JOB_NAME  -lMH_min 13.5 -lMH_max 13.7 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin7'$SLURM_JOB_NAME  -lMH_min 13.7 -lMH_max 13.9 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin8'$SLURM_JOB_NAME  -lMH_min 13.9 -lMH_max 14.1 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin9'$SLURM_JOB_NAME  -lMH_min 14.1 -lMH_max 14.3 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
srun python -u forGroup_profile_CS82.py -sample 'Mbin10'$SLURM_JOB_NAME -lMH_min 14.3 -lMH_max 15.5 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 50
                                          
