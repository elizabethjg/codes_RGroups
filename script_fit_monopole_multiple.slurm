#!/bin/bash
#SBATCH --job-name=fitMbins
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
##SBATCH --tasks-per-node=20
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/

export OMP_NUM_THREADS=56
export MKL_NUM_THREADS=56

source activate py2env

case $SLURM_NODEID in
    0) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin1_2.fits' -ncores 10;;
    1) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin2_2.fits' -ncores 10;;
    2) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin3_2.fits' -ncores 10;;
    3) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin4_2.fits' -ncores 10;;
    4) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin5_2.fits' -ncores 10;;
    5) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin6_2.fits' -ncores 10;;
    6) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin7_2.fits' -ncores 10;;
    7) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin8_2.fits' -ncores 10;;
    8) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin9_2.fits' -ncores 10;;
    9) srun -u python fit_profile_monopole_misscentred_pcc.py \
    -folder '/mnt/clemente/lensing/RodriguezGroups/N_all/' \
    -file 'profile_Mbin10_2.fits' -ncores 10;;
esac 
