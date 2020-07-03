#!/bin/bash
#SBATCH --job-name=Mbin_cM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=12:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar

. /etc/profile
source $HOME/.bashrc

cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

conda activate py2env

srun python -u forGroup_profile.py -sample 'Mbin1_cM'  -lMH_min 12.5 -lMH_max 12.7 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin2_cM'  -lMH_min 12.7 -lMH_max 12.9 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin3_cM'  -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin4_cM'  -lMH_min 13.1 -lMH_max 13.3 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin5_cM'  -lMH_min 13.3 -lMH_max 13.5 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin6_cM'  -lMH_min 13.5 -lMH_max 13.7 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin7_cM'  -lMH_min 13.7 -lMH_max 13.9 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin8_cM'  -lMH_min 13.9 -lMH_max 14.1 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin9_cM'  -lMH_min 14.1 -lMH_max 14.3 -C_BG_min 0.3 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin10_cM' -lMH_min 14.3 -lMH_max 15.5 -C_BG_min 0.3 -z_max 0.2 -ncores 56
                                                                                       
srun python -u forGroup_profile.py -sample 'Mbin1_zL_cM'  -lMH_min 12.5 -lMH_max 12.7 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin2_zL_cM'  -lMH_min 12.7 -lMH_max 12.9 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin3_zL_cM'  -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin4_zL_cM'  -lMH_min 13.1 -lMH_max 13.3 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin5_zL_cM'  -lMH_min 13.3 -lMH_max 13.5 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin6_zL_cM'  -lMH_min 13.5 -lMH_max 13.7 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin7_zL_cM'  -lMH_min 13.7 -lMH_max 13.9 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin8_zL_cM'  -lMH_min 13.9 -lMH_max 14.1 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin9_zL_cM'  -lMH_min 14.1 -lMH_max 14.3 -C_BG_min 0.3 -z_max 0.13 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin10_zL_cM' -lMH_min 14.3 -lMH_max 15.5 -C_BG_min 0.3 -z_max 0.13 -ncores 56
                                                                                       
srun python -u forGroup_profile.py -sample 'Mbin1_zH_cM'  -lMH_min 12.5 -lMH_max 12.7 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin2_zH_cM'  -lMH_min 12.7 -lMH_max 12.9 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin3_zH_cM'  -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin4_zH_cM'  -lMH_min 13.1 -lMH_max 13.3 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin5_zH_cM'  -lMH_min 13.3 -lMH_max 13.5 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin6_zH_cM'  -lMH_min 13.5 -lMH_max 13.7 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin7_zH_cM'  -lMH_min 13.7 -lMH_max 13.9 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin8_zH_cM'  -lMH_min 13.9 -lMH_max 14.1 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin9_zH_cM'  -lMH_min 14.1 -lMH_max 14.3 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin10_zH_cM' -lMH_min 14.3 -lMH_max 15.5 -C_BG_min 0.3 -z_min 0.13 -z_max 0.2 -ncores 56
