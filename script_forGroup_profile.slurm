#!/bin/bash
#SBATCH --job-name=Mbin
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

srun python -u forGroup_profile.py -sample 'Mbin1_cM_mzL'  -lMH_min 12.5 -lMH_max 12.7 -z_max 0.134 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin2_cM_mzL'  -lMH_min 12.7 -lMH_max 12.9 -z_max 0.149 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin3_cM_mzL'  -lMH_min 12.9 -lMH_max 13.1 -z_max 0.163 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin4_cM_mzL'  -lMH_min 13.1 -lMH_max 13.3 -z_max 0.166 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin5_cM_mzL'  -lMH_min 13.3 -lMH_max 13.5 -z_max 0.166 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin6_cM_mzL'  -lMH_min 13.5 -lMH_max 13.7 -z_max 0.159 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin7_cM_mzL'  -lMH_min 13.7 -lMH_max 13.9 -z_max 0.133 -C_BG_min 2.73 -C_BG_max 100. -ncores 40
srun python -u forGroup_profile.py -sample 'Mbin8_cM_mzL'  -lMH_min 13.9 -lMH_max 14.1 -z_max 0.125 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin9_cM_mzL'  -lMH_min 14.1 -lMH_max 14.3 -z_max 0.116 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin10_cM_mzL' -lMH_min 14.3 -lMH_max 15.5 -z_max 0.111 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
                                                                                                
srun python -u forGroup_profile.py -sample 'Mbin1_cM_mzH'  -lMH_min 12.5 -lMH_max 12.7 -z_min 0.134 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin2_cM_mzH'  -lMH_min 12.7 -lMH_max 12.9 -z_min 0.149 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin3_cM_mzH'  -lMH_min 12.9 -lMH_max 13.1 -z_min 0.163 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin4_cM_mzH'  -lMH_min 13.1 -lMH_max 13.3 -z_min 0.166 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin5_cM_mzH'  -lMH_min 13.3 -lMH_max 13.5 -z_min 0.166 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin6_cM_mzH'  -lMH_min 13.5 -lMH_max 13.7 -z_min 0.159 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin7_cM_mzH'  -lMH_min 13.7 -lMH_max 13.9 -z_min 0.133 -C_BG_min 2.73 -C_BG_max 100. -ncores 40
srun python -u forGroup_profile.py -sample 'Mbin8_cM_mzH'  -lMH_min 13.9 -lMH_max 14.1 -z_min 0.125 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin9_cM_mzH'  -lMH_min 14.1 -lMH_max 14.3 -z_min 0.116 -C_BG_min 2.73 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin10_cM_mzH' -lMH_min 14.3 -lMH_max 15.5 -z_min 0.111 -C_BG_min 2.73 -C_BG_max 100. -ncores 46

srun python -u forGroup_profile.py -sample 'Mbin1_mzL'  -lMH_min 12.5 -lMH_max 12.7 -z_max 0.135 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin2_mzL'  -lMH_min 12.7 -lMH_max 12.9 -z_max 0.150 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin3_mzL'  -lMH_min 12.9 -lMH_max 13.1 -z_max 0.161 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin4_mzL'  -lMH_min 13.1 -lMH_max 13.3 -z_max 0.164 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin5_mzL'  -lMH_min 13.3 -lMH_max 13.5 -z_max 0.163 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin6_mzL'  -lMH_min 13.5 -lMH_max 13.7 -z_max 0.156 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin7_mzL'  -lMH_min 13.7 -lMH_max 13.9 -z_max 0.134 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin8_mzL'  -lMH_min 13.9 -lMH_max 14.1 -z_max 0.126 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin9_mzL'  -lMH_min 14.1 -lMH_max 14.3 -z_max 0.118 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin10_mzL' -lMH_min 14.3 -lMH_max 15.5 -z_max 0.111 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
                                                                                                                                      
srun python -u forGroup_profile.py -sample 'Mbin1_mzH'  -lMH_min 12.5 -lMH_max 12.7 -z_min 0.135 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin2_mzH'  -lMH_min 12.7 -lMH_max 12.9 -z_min 0.150 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin3_mzH'  -lMH_min 12.9 -lMH_max 13.1 -z_min 0.161 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin4_mzH'  -lMH_min 13.1 -lMH_max 13.3 -z_min 0.164 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin5_mzH'  -lMH_min 13.3 -lMH_max 13.5 -z_min 0.163 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin6_mzH'  -lMH_min 13.5 -lMH_max 13.7 -z_min 0.156 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin7_mzH'  -lMH_min 13.7 -lMH_max 13.9 -z_min 0.134 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin8_mzH'  -lMH_min 13.9 -lMH_max 14.1 -z_min 0.126 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin9_mzH'  -lMH_min 14.1 -lMH_max 14.3 -z_min 0.118 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
srun python -u forGroup_profile.py -sample 'Mbin10_mzH' -lMH_min 14.3 -lMH_max 15.5 -z_min 0.111 -C_BG_min 0.0 -C_BG_max 100. -ncores 46
