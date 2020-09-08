#!/bin/bash
#SBATCH --job-name=Nbin1_z
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
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

##srun python -u forGroup_profile.py -sample 'N1_Mbin1'   -N_min 1 -N_max 2   -lMH_min 12.5 -lMH_max 12.9 -C_BG_min 0. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N1_Mbin2'   -N_min 1 -N_max 2   -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 0. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N1_Mbin3'   -N_min 1 -N_max 2   -lMH_min 13.1 -lMH_max 13.5 -C_BG_min 0. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N23_Mbin1'  -N_min 2 -N_max 4   -lMH_min 12.5 -lMH_max 13.5 -C_BG_min 0. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N23_Mbin2'  -N_min 2 -N_max 4   -lMH_min 13.5 -lMH_max 14.5 -C_BG_min 0. -z_max 0.2 -ncores 32
srun python -u forGroup_profile.py -sample 'N4M_Mbin1'  -N_min 4 -N_max 100 -lMH_min 12.5 -lMH_max 13.7 -C_BG_min 0. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'N4M_Mbin2'  -N_min 4 -N_max 100 -lMH_min 13.7 -lMH_max 14.2 -C_BG_min 0. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'N4M_Mbin3'  -N_min 4 -N_max 100 -lMH_min 14.2 -lMH_max 15.5 -C_BG_min 0. -z_max 0.2 -ncores 56
  
##srun python -u forGroup_profile.py -sample 'N1_Mbin1_cM'   -N_min 1 -N_max 2   -lMH_min 12.5 -lMH_max 12.9 -C_BG_min 3. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N1_Mbin2_cM'   -N_min 1 -N_max 2   -lMH_min 12.9 -lMH_max 13.1 -C_BG_min 3. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N1_Mbin3_cM'   -N_min 1 -N_max 2   -lMH_min 13.1 -lMH_max 13.5 -C_BG_min 3. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N23_Mbin1_cM'  -N_min 2 -N_max 4   -lMH_min 12.5 -lMH_max 13.5 -C_BG_min 3. -z_max 0.2 -ncores 32
##srun python -u forGroup_profile.py -sample 'N23_Mbin2_cM'  -N_min 2 -N_max 4   -lMH_min 13.5 -lMH_max 14.5 -C_BG_min 3. -z_max 0.2 -ncores 32
srun python -u forGroup_profile.py -sample 'N4M_Mbin1_cM'  -N_min 4 -N_max 100 -lMH_min 12.5 -lMH_max 13.7 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'N4M_Mbin2_cM'  -N_min 4 -N_max 100 -lMH_min 13.7 -lMH_max 14.2 -C_BG_min 3. -z_max 0.2 -ncores 56
srun python -u forGroup_profile.py -sample 'N4M_Mbin3_cM'  -N_min 4 -N_max 100 -lMH_min 14.2 -lMH_max 15.5 -C_BG_min 3. -z_max 0.2 -ncores 56

