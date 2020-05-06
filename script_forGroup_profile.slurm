#!/bin/bash
#SBATCH --job-name=Mbins
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

srun python -u forGroup_profile.py -sample 'Mbin1' -lMH_min 12. -lMH_max 12.4 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin2' -lMH_min 12.4 -lMH_max 12.8 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin3' -lMH_min 12.8 -lMH_max 13.2 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin4' -lMH_min 13.2 -lMH_max 13.6 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin5' -lMH_min 13.6 -lMH_max 14.0 -ncores 56
srun python -u forGroup_profile.py -sample 'Mbin6' -lMH_min 14.0 -lMH_max 15.5 -ncores 56



##srun python -u forGroup_profile.py -sample 'Mbin1' -lMH_min 0. -lMH_max 12.0310770620 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin2' -lMH_min 12.0310770620 -lMH_max 12.2896221699 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin3' -lMH_min 12.2896221699 -lMH_max 12.5548037732 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin4' -lMH_min 12.5548037732 -lMH_max 12.9098106480 -ncores 56
##srun python -u forGroup_profile.py -sample 'Mbin5' -lMH_min 12.9098106480 -lMH_max 15.5 -ncores 56

##srun python -u forGroup_profile.py -sample 'N1_conM' -N_min 1 -N_max 2 -C_BG_min 2.73 -ncores 56
 
