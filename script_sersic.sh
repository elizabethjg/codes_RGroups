#!/bin/bash
cd /mnt/clemente/lensing/RodriguezGroups/codes_RGroups/
conda activate py2env

python -u forGroup_profile.py -sample 'Mbin1_n1_cM' -lMH_min 12.5 -lMH_max 12.7 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48
python -u forGroup_profile.py -sample 'Mbin2_n1_cM' -lMH_min 12.7 -lMH_max 12.9 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48
python -u forGroup_profile.py -sample 'Mbin3_n1_cM' -lMH_min 12.9 -lMH_max 13.1 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48
python -u forGroup_profile.py -sample 'Mbin4_n1_cM' -lMH_min 13.1 -lMH_max 13.3 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48
python -u forGroup_profile.py -sample 'Mbin5_n1_cM' -lMH_min 13.3 -lMH_max 13.5 -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48
python -u forGroup_profile.py -sample 'Mbin6_n1_cM' -lMH_min 13.5 -lMH_max 15.0   -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48

##python -u forGroup_profile.py -sample 'Mbin1_n1_mzL' -lMH_min 12.5 -lMH_max 12.7 -z_max 0.141 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n1_mzL' -lMH_min 12.7 -lMH_max 12.9 -z_max 0.166 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n1_mzL' -lMH_min 12.9 -lMH_max 13.1 -z_max 0.187 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin4_n1_mzL' -lMH_min 13.1 -lMH_max 13.3 -z_max 0.210 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin5_n1_mzL' -lMH_min 13.3 -lMH_max 13.5 -z_max 0.242 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin6_n1_mzL' -lMH_min 13.5 -lMH_max 15.0   -z_max 0.259 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##
##                                       
##python -u forGroup_profile.py -sample 'Mbin1_n1_mzL_cM' -lMH_min 12.5 -lMH_max 12.7 -z_max 0.141 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n1_mzL_cM' -lMH_min 12.7 -lMH_max 12.9 -z_max 0.166 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n1_mzL_cM' -lMH_min 12.9 -lMH_max 13.1 -z_max 0.187 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin4_n1_mzL_cM' -lMH_min 13.1 -lMH_max 13.3 -z_max 0.210 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin5_n1_mzL_cM' -lMH_min 13.3 -lMH_max 13.5 -z_max 0.242 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin6_n1_mzL_cM' -lMH_min 13.5 -lMH_max 15.0   -z_max 0.259 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##
##                                                                                                                                                   
##python -u forGroup_profile.py -sample 'Mbin1_n1_mzH_cM' -lMH_min 12.5 -lMH_max 12.7 -z_min 0.141 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n1_mzH_cM' -lMH_min 12.7 -lMH_max 12.9 -z_min 0.166 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n1_mzH_cM' -lMH_min 12.9 -lMH_max 13.1 -z_min 0.187 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin4_n1_mzH_cM' -lMH_min 13.1 -lMH_max 13.3 -z_min 0.210 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin5_n1_mzH_cM' -lMH_min 13.3 -lMH_max 13.5 -z_min 0.242 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin6_n1_mzH_cM' -lMH_min 13.5 -lMH_max 15.0   -z_min 0.259 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                       
##python -u forGroup_profile.py -sample 'Mbin1_n23_cM' -lMH_min 12.5 -lMH_max 13.1 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n23_cM' -lMH_min 13.1 -lMH_max 13.7 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n23_cM' -lMH_min 13.7 -lMH_max 14.5 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'Mbin1_n23_mzL_cM' -lMH_min 12.5 -lMH_max 13.1 -z_max 0.099 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n23_mzL_cM' -lMH_min 13.1 -lMH_max 13.7 -z_max 0.140 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n23_mzL_cM' -lMH_min 13.7 -lMH_max 14.5 -z_max 0.184 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                                                                                                                                            
##python -u forGroup_profile.py -sample 'Mbin1_n23_mzH_cM' -lMH_min 12.5 -lMH_max 13.1 -z_min 0.099 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n23_mzH_cM' -lMH_min 13.1 -lMH_max 13.7 -z_min 0.140 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n23_mzH_cM' -lMH_min 13.7 -lMH_max 14.5 -z_min 0.184 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'Mbin1_n4M_cM' -lMH_min 12.5 -lMH_max 13.3 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n4M_cM' -lMH_min 13.3 -lMH_max 13.7 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n4M_cM' -lMH_min 13.7 -lMH_max 15.5 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                         
##python -u forGroup_profile.py -sample 'Mbin1_n4M_mzL_cM' -lMH_min 12.5 -lMH_max 13.3 -z_max 0.076 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n4M_mzL_cM' -lMH_min 13.3 -lMH_max 13.7 -z_max 0.093 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n4M_mzL_cM' -lMH_min 13.7 -lMH_max 15.5 -z_max 0.117 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                                                                                                                   
##python -u forGroup_profile.py -sample 'Mbin1_n4M_mzH_cM' -lMH_min 12.5 -lMH_max 13.3 -z_min 0.076 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin2_n4M_mzH_cM' -lMH_min 13.3 -lMH_max 13.7 -z_min 0.093 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'Mbin3_n4M_mzH_cM' -lMH_min 13.7 -lMH_max 15.5 -z_min 0.117 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'N1_Mbin1'  -N_max 2 -lMH_min 12.5 -lMH_max 12.7 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin2'  -N_max 2 -lMH_min 12.7 -lMH_max 12.9 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin3'  -N_max 2 -lMH_min 12.9 -lMH_max 13.1 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin4'  -N_max 2 -lMH_min 13.1 -lMH_max 13.3 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin5'  -N_max 2 -lMH_min 13.3 -lMH_max 13.5 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin6'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0 -z_max 1. -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin6_cM'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_max 1. -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                     
##python -u forGroup_profile.py -sample 'N1_Mbin1_mzL'  -N_max 2 -lMH_min 12.5 -lMH_max 12.7 -z_max 0.141 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin2_mzL'  -N_max 2 -lMH_min 12.7 -lMH_max 12.9 -z_max 0.166 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin3_mzL'  -N_max 2 -lMH_min 12.9 -lMH_max 13.1 -z_max 0.187 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin4_mzL'  -N_max 2 -lMH_min 13.1 -lMH_max 13.3 -z_max 0.210 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin5_mzL'  -N_max 2 -lMH_min 13.3 -lMH_max 13.5 -z_max 0.242 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin6_mzL'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_max 0.259 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin6_cM_mzL'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_max 0.259 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'N1_Mbin1_mzH'  -N_max 2 -lMH_min 12.5 -lMH_max 12.7 -z_min 0.141 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin2_mzH'  -N_max 2 -lMH_min 12.7 -lMH_max 12.9 -z_min 0.166 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin3_mzH'  -N_max 2 -lMH_min 12.9 -lMH_max 13.1 -z_min 0.187 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin4_mzH'  -N_max 2 -lMH_min 13.1 -lMH_max 13.3 -z_min 0.210 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin5_mzH'  -N_max 2 -lMH_min 13.3 -lMH_max 13.5 -z_min 0.242 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin6_mzH'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_min 0.259 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N1_Mbin6_cM_mzH'  -N_max 2 -lMH_min 13.5 -lMH_max 15.0   -z_min 0.259 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                     
##python -u forGroup_profile.py -sample 'N2-3_Mbin1'  -N_min 2 -N_max 4 -lMH_min 12.5 -lMH_max 13.1 -z_max 1.0 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N2-3_Mbin2'  -N_min 2 -N_max 4 -lMH_min 13.1 -lMH_max 13.7 -z_max 1.0 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N2-3_Mbin3'  -N_min 2 -N_max 4 -lMH_min 13.7 -lMH_max 14.5 -z_max 1.0 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'N2-3_Mbin1_mzL'  -N_min 2 -N_max 4 -lMH_min 12.5 -lMH_max 13.1 -z_max 0.099 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N2-3_Mbin2_mzL'  -N_min 2 -N_max 4 -lMH_min 13.1 -lMH_max 13.7 -z_max 0.140 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N2-3_Mbin3_mzL'  -N_min 2 -N_max 4 -lMH_min 13.7 -lMH_max 14.5 -z_max 0.184 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##                                                                                                                                                            
##python -u forGroup_profile.py -sample 'N2-3_Mbin1_mzH'  -N_min 2 -N_max 4 -lMH_min 12.5 -lMH_max 13.1 -z_min 0.099 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N2-3_Mbin2_mzH'  -N_min 2 -N_max 4 -lMH_min 13.1 -lMH_max 13.7 -z_min 0.140 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N2-3_Mbin3_mzH'  -N_min 2 -N_max 4 -lMH_min 13.7 -lMH_max 14.5 -z_min 0.184 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'N4M_Mbin1'  -N_min 4 -lMH_min 12.5 -lMH_max 13.3 -z_max 1.0 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin2'  -N_min 4 -lMH_min 13.3 -lMH_max 13.7 -z_max 1.0 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin3'  -N_min 4 -lMH_min 13.7 -lMH_max 15.5 -z_max 1.0 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##                                         
##python -u forGroup_profile.py -sample 'N4M_Mbin1_mzL'  -N_min 4 -lMH_min 12.5 -lMH_max 13.3 -z_max 0.076 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin2_mzL'  -N_min 4 -lMH_min 13.3 -lMH_max 13.7 -z_max 0.093 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin3_mzL'  -N_min 4 -lMH_min 13.7 -lMH_max 15.5 -z_max 0.117 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##                                                                                                                                                           
##python -u forGroup_profile.py -sample 'N4M_Mbin1_mzH'  -N_min 4 -lMH_min 12.5 -lMH_max 13.3 -z_min 0.076 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin2_mzH'  -N_min 4 -lMH_min 13.3 -lMH_max 13.7 -z_min 0.093 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin3_mzH'  -N_min 4 -lMH_min 13.7 -lMH_max 15.5 -z_min 0.117 -C_BG_min 0.0 -C_BG_max 100. -ncores 48
##
##python -u forGroup_profile.py -sample 'N4M_Mbin1_cM'  -N_min 4 -lMH_min 12.5 -lMH_max 13.3 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin2_cM'  -N_min 4 -lMH_min 13.3 -lMH_max 13.7 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin3_cM'  -N_min 4 -lMH_min 13.7 -lMH_max 15.5 -z_max 1.0 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                         
##python -u forGroup_profile.py -sample 'N4M_Mbin1_mzL_cM'  -N_min 4 -lMH_min 12.5 -lMH_max 13.3 -z_max 0.076 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin2_mzL_cM'  -N_min 4 -lMH_min 13.3 -lMH_max 13.7 -z_max 0.093 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin3_mzL_cM'  -N_min 4 -lMH_min 13.7 -lMH_max 15.5 -z_max 0.117 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##                                                                                                                                                           
##python -u forGroup_profile.py -sample 'N4M_Mbin1_mzH_cM'  -N_min 4 -lMH_min 12.5 -lMH_max 13.3 -z_min 0.076 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin2_mzH_cM'  -N_min 4 -lMH_min 13.3 -lMH_max 13.7 -z_min 0.093 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
##python -u forGroup_profile.py -sample 'N4M_Mbin3_mzH_cM'  -N_min 4 -lMH_min 13.7 -lMH_max 15.5 -z_min 0.117 -C_BG_min 2.73 -C_BG_max 100. -ncores 48
