#!/bin/bash
#SBATCH --job-name=Mbin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=46
#SBATCH --time=6:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=46
export MKL_NUM_THREADS=46

source activate py2env

srun python -u forGroup_profile.py -sample 'N1_Mbin1_cM'  -N_max 2 -lMH_min 12.5 -lMH_max 12.7 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin2_cM'  -N_max 2 -lMH_min 12.7 -lMH_max 12.9 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin3_cM'  -N_max 2 -lMH_min 12.9 -lMH_max 13.1 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin4_cM'  -N_max 2 -lMH_min 13.1 -lMH_max 13.3 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin5_cM'  -N_max 2 -lMH_min 13.3 -lMH_max 13.5 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 46

srun python -u forGroup_profile.py -sample 'N1_Mbin1_cM_zL'  -N_max 2 -lMH_min 12.5 -lMH_max 12.7 -z_max 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin2_cM_zL'  -N_max 2 -lMH_min 12.7 -lMH_max 12.9 -z_max 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin3_cM_zL'  -N_max 2 -lMH_min 12.9 -lMH_max 13.1 -z_max 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin4_cM_zL'  -N_max 2 -lMH_min 13.1 -lMH_max 13.3 -z_max 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin5_cM_zL'  -N_max 2 -lMH_min 13.3 -lMH_max 13.5 -z_max 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                                            
srun python -u forGroup_profile.py -sample 'N1_Mbin1_cM_zH'  -N_max 2 -lMH_min 12.5 -lMH_max 12.7 -z_min 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin2_cM_zH'  -N_max 2 -lMH_min 12.7 -lMH_max 12.9 -z_min 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin3_cM_zH'  -N_max 2 -lMH_min 12.9 -lMH_max 13.1 -z_min 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin4_cM_zH'  -N_max 2 -lMH_min 13.1 -lMH_max 13.3 -z_min 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N1_Mbin5_cM_zH'  -N_max 2 -lMH_min 13.3 -lMH_max 13.5 -z_min 0.167 -C_BG_min 2.73 -C_BG_max 100. -ncores 46

srun python -u forGroup_profile.py -sample 'N2-3_Mbin1_cM'  -N_min 2 -N_max 4 -lMH_min 12.5 -lMH_max 12.9 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin2_cM'  -N_min 2 -N_max 4 -lMH_min 12.9 -lMH_max 13.3 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin3_cM'  -N_min 2 -N_max 4 -lMH_min 13.3 -lMH_max 13.7 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin4_cM'  -N_min 2 -N_max 4 -lMH_min 13.7 -lMH_max 14.5 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                      
srun python -u forGroup_profile.py -sample 'N2-3_Mbin1_cM_zL'  -N_min 2 -N_max 4 -lMH_min 12.5 -lMH_max 12.9 -z_max 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin2_cM_zL'  -N_min 2 -N_max 4 -lMH_min 12.9 -lMH_max 13.3 -z_max 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin3_cM_zL'  -N_min 2 -N_max 4 -lMH_min 13.3 -lMH_max 13.7 -z_max 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin4_cM_zL'  -N_min 2 -N_max 4 -lMH_min 13.7 -lMH_max 14.5 -z_max 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                               
srun python -u forGroup_profile.py -sample 'N2-3_Mbin1_cM_zH'  -N_min 2 -N_max 4 -lMH_min 12.5 -lMH_max 12.9 -z_min 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin2_cM_zH'  -N_min 2 -N_max 4 -lMH_min 12.9 -lMH_max 13.3 -z_min 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin3_cM_zH'  -N_min 2 -N_max 4 -lMH_min 13.3 -lMH_max 13.7 -z_min 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N2-3_Mbin4_cM_zH'  -N_min 2 -N_max 4 -lMH_min 13.7 -lMH_max 14.5 -z_min 0.118 -C_BG_min 2.73 -C_BG_max 100. -ncores 46

srun python -u forGroup_profile.py -sample 'N4-10_Mbin1_cM'  -N_min 4 -N_max 11 -lMH_min 12.5 -lMH_max 13.0 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N4-10_Mbin2_cM'  -N_min 4 -N_max 11 -lMH_min 13.0 -lMH_max 13.5 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N4-10_Mbin3_cM'  -N_min 4 -N_max 11 -lMH_min 13.5 -lMH_max 14.5 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                      
srun python -u forGroup_profile.py -sample 'N4-10_Mbin1_cM_zL'  -N_min 4 -N_max 11 -lMH_min 12.5 -lMH_max 13.0 -z_max 0.096 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N4-10_Mbin2_cM_zL'  -N_min 4 -N_max 11 -lMH_min 13.0 -lMH_max 13.5 -z_max 0.096 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N4-10_Mbin3_cM_zL'  -N_min 4 -N_max 11 -lMH_min 13.5 -lMH_max 14.5 -z_max 0.096 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                                                                        
srun python -u forGroup_profile.py -sample 'N4-10_Mbin1_cM_zH'  -N_min 4 -N_max 11 -lMH_min 12.5 -lMH_max 13.0 -z_min 0.096 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N4-10_Mbin2_cM_zH'  -N_min 4 -N_max 11 -lMH_min 13.0 -lMH_max 13.5 -z_min 0.096 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N4-10_Mbin3_cM_zH'  -N_min 4 -N_max 11 -lMH_min 13.5 -lMH_max 14.5 -z_min 0.096 -C_BG_min 2.73 -C_BG_max 100. -ncores 46

srun python -u forGroup_profile.py -sample 'N10M_Mbin1_cM'  -N_min 10 -lMH_min 13.0 -lMH_max 14.0 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N10M_Mbin2_cM'  -N_min 10 -lMH_min 14.0 -lMH_max 15.0 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                      
srun python -u forGroup_profile.py -sample 'N10M_Mbin1_cM_zL'  -N_min 10 -lMH_min 13.0 -lMH_max 14.0 -z_max 0.083 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N10M_Mbin2_cM_zL'  -N_min 10 -lMH_min 14.0 -lMH_max 15.0 -z_max 0.083 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                                                               
srun python -u forGroup_profile.py -sample 'N10M_Mbin1_cM_zH'  -N_min 10 -lMH_min 13.0 -lMH_max 14.0 -z_min 0.083 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'N10M_Mbin2_cM_zH'  -N_min 10 -lMH_min 14.0 -lMH_max 15.0 -z_min 0.083 -C_BG_min 2.73 -C_BG_max 100. -ncores 46

##srun python -u forGroup_profile.py -sample 'Mbin1_zL'  -lMH_min 12.5 -lMH_max 12.7 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin2_zL'  -lMH_min 12.7 -lMH_max 12.9 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin3_zL'  -lMH_min 12.9 -lMH_max 13.1 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin4_zL'  -lMH_min 13.1 -lMH_max 13.3 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin5_zL'  -lMH_min 13.3 -lMH_max 13.5 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin6_zL'  -lMH_min 13.5 -lMH_max 13.7 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin7_zL'  -lMH_min 13.7 -lMH_max 13.9 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin8_zL'  -lMH_min 13.9 -lMH_max 14.1 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin9_zL'  -lMH_min 14.1 -lMH_max 14.3 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin10_zL' -lMH_min 14.3 -lMH_max 15.5 -z_max 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
                                                                                                                  
##srun python -u forGroup_profile.py -sample 'Mbin1_zH'  -lMH_min 12.5 -lMH_max 12.7 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin2_zH'  -lMH_min 12.7 -lMH_max 12.9 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin3_zH'  -lMH_min 12.9 -lMH_max 13.1 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin4_zH'  -lMH_min 13.1 -lMH_max 13.3 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin5_zH'  -lMH_min 13.3 -lMH_max 13.5 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin6_zH'  -lMH_min 13.5 -lMH_max 13.7 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin7_zH'  -lMH_min 13.7 -lMH_max 13.9 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin8_zH'  -lMH_min 13.9 -lMH_max 14.1 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin9_zH'  -lMH_min 14.1 -lMH_max 14.3 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
##srun python -u forGroup_profile.py -sample 'Mbin10_zH' -lMH_min 14.3 -lMH_max 15.5 -z_min 0.15 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
