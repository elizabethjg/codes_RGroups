#!/bin/bash
#SBATCH --job-name=Mbin_cM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=6:00:00
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

srun python -u forGroup_profile.py -sample 'Mbin1_cM'  -lMH_min 12.5 -lMH_max 13.0 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin2_cM'  -lMH_min 13.0 -lMH_max 13.3 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin3_cM'  -lMH_min 13.3 -lMH_max 13.6 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin4_cM'  -lMH_min 13.6 -lMH_max 13.9 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin5_cM'  -lMH_min 13.9 -lMH_max 14.2 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin6_cM'  -lMH_min 14.2 -lMH_max 15.5 -C_BG_min 3. -z_max 0.2 -ncores 56

##srun python -u forGroup_profile.py -sample 'Mbin1_zL'  -lMH_min 12.5 -lMH_max 12.9 -C_BG_min 0. -z_max 0.13 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin2_zL'  -lMH_min 12.9 -lMH_max 13.3 -C_BG_min 0. -z_max 0.13 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin3_zL'  -lMH_min 13.3 -lMH_max 13.7 -C_BG_min 0. -z_max 0.13 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin4_zL'  -lMH_min 13.7 -lMH_max 14.1 -C_BG_min 0. -z_max 0.13 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin5_zL'  -lMH_min 14.1 -lMH_max 14.5 -C_BG_min 0. -z_max 0.13 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin6_zL'  -lMH_min 14.5 -lMH_max 15.5 -C_BG_min 0. -z_max 0.13 -ncores 45
                                                                                                  
##srun python -u forGroup_profile.py -sample 'Mbin1_zH'  -lMH_min 12.5 -lMH_max 12.9 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin2_zH'  -lMH_min 12.9 -lMH_max 13.3 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin3_zH'  -lMH_min 13.3 -lMH_max 13.7 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin4_zH'  -lMH_min 13.7 -lMH_max 14.1 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin5_zH'  -lMH_min 14.1 -lMH_max 14.5 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 45
##srun python -u forGroup_profile.py -sample 'Mbin6_zH'  -lMH_min 14.5 -lMH_max 15.5 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 45


                                                                                       
##srun python -u forGroup_profile.py -sample 'Mbin1_zL'  -lMH_min 12.5 -lMH_max 12.7 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin2_zL'  -lMH_min 12.7 -lMH_max 12.9 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin3_zL'  -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin4_zL'  -lMH_min 13.1 -lMH_max 13.3 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin5_zL'  -lMH_min 13.3 -lMH_max 13.5 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin6_zL'  -lMH_min 13.5 -lMH_max 13.7 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin7_zL'  -lMH_min 13.7 -lMH_max 13.9 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin8_zL'  -lMH_min 13.9 -lMH_max 14.1 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin9_zL'  -lMH_min 14.1 -lMH_max 14.3 -C_BG_min 0. -z_max 0.13 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin10_zL' -lMH_min 14.3 -lMH_max 15.5 -C_BG_min 0. -z_max 0.13 -ncores 50
##                                                                                             
##srun python -u forGroup_profile.py -sample 'Mbin1_zH'  -lMH_min 12.5 -lMH_max 12.7 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin2_zH'  -lMH_min 12.7 -lMH_max 12.9 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin3_zH'  -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin4_zH'  -lMH_min 13.1 -lMH_max 13.3 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin5_zH'  -lMH_min 13.3 -lMH_max 13.5 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin6_zH'  -lMH_min 13.5 -lMH_max 13.7 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin7_zH'  -lMH_min 13.7 -lMH_max 13.9 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin8_zH'  -lMH_min 13.9 -lMH_max 14.1 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin9_zH'  -lMH_min 14.1 -lMH_max 14.3 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
##srun python -u forGroup_profile.py -sample 'Mbin10_zH' -lMH_min 14.3 -lMH_max 15.5 -C_BG_min 0. -z_min 0.13 -z_max 0.2 -ncores 50
