#!/bin/bash
#SBATCH --job-name=FOF_redo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=56
export MKL_NUM_THREADS=56

source activate py2env

srun python -u forGroup_profile.py -sample 'Mbin7' -lMH_min 13.7 -lMH_max 13.9  -C_BG_min 0. -C_BG_max 100. -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin10' -lMH_min 14.3 -lMH_max 15.5  -C_BG_min 0. -C_BG_max 100. -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin7_conM' -lMH_min 13.7 -lMH_max 13.9  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin10_conM' -lMH_min 14.3 -lMH_max 15.5  -C_BG_min 2.73 -C_BG_max 100. -ncores 56

##srun python -u forGroup_profile.py -sample 'Mbin1_conM' -lMH_min 12.5 -lMH_max 12.7  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin2_conM' -lMH_min 12.7 -lMH_max 12.9  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin3_conM' -lMH_min 12.9 -lMH_max 13.1  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin4_conM' -lMH_min 13.1 -lMH_max 13.3  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin5_conM' -lMH_min 13.3 -lMH_max 13.5  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin6_conM' -lMH_min 13.5 -lMH_max 13.7  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin7_conM' -lMH_min 13.9 -lMH_max 14.1  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin8_conM' -lMH_min 14.1 -lMH_max 14.3  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin9_conM' -lMH_min 14.3 -lMH_max 14.5  -C_BG_min 2.73 -C_BG_max 100. -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin10_conM' -lMH_min 14.5 -lMH_max 15.5 -C_BG_min 2.73 -C_BG_max 100. -ncores 56



##srun python -u forGroup_profile.py -sample 'Mbin1' -lMH_min 0. -lMH_max 12.0310770620 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin2' -lMH_min 12.0310770620 -lMH_max 12.2896221699 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin3' -lMH_min 12.2896221699 -lMH_max 12.5548037732 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin4' -lMH_min 12.5548037732 -lMH_max 12.9098106480 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin5' -lMH_min 12.9098106480 -lMH_max 15.5 -ncores 56

##srun python -u forGroup_profile.py -sample 'N1_conM' -N_min 1 -N_max 2 -C_BG_min 2.73 -ncores 56
 
