#!/bin/bash
#SBATCH --job-name=fitMbins
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --tasks-per-node=20
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=20
export MKL_NUM_THREADS=20

source activate py2env

case $SLURM_NODEID in
    0) srun -u python forGroup_profile.py -sample 'pru1' \
    -lMH_min 0. -lMH_max 11.8 -ncores 10 ;;
    1) srun -u python forGroup_profile.py -sample 'pru2' \
    -lMH_min 12. -lMH_max 12.5 -ncores 10 ;;
