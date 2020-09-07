import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/eli/python_codes')
import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits
import corner
from profiles_fit import *
import os
import argparse



def make_plot_misscentred_monopole(file_name,folder,ax,lab):
          
     
     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data
     
     file_mcmc = 'monopole_pcconly_'+file_name[:-4]+'out'
          
     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Rmean   = h['RADIUS_HALO_mean']
     ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))/0.7
     soff = 0.4*Rmean
     Vdisp   = h['VDISP_DYN_mean']
     
     mcmc = (np.loadtxt(folder+file_mcmc)).T
          
     labels = ['M200','pcc']
     
     mout    = np.percentile(mcmc[0][1000:], [16, 50, 84])
     pcc_out = np.percentile(mcmc[1][1000:], [16, 50, 84])
     
     
          
     # e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     
     cmin, cmax = h['C_BG_min'],h['C_BG_max']
     Mmin, Mmax = h['lMH_min'],h['lMH_max']
     zmin, zmax = h['z_min'],h['z_max']
     
     print 'C_BG, min, max'
     print cmin, cmax
     print 'lMH, min, max'
     print Mmin, Mmax
     print 'z, min, max'
     print zmin, zmax
     
     
          
     M200   = 10**(mout[1])
     e_M200 = (10**(mout[1])*np.log(10.)*np.diff(mout))
     
     
     #--------------------
          
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
          
     multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                                   misscentred = True,s_off = soff,
                                   ellip=0,z=zmean,components = ['t'],
                                   verbose=False,ncores=2)
          
     Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])     
          
     Gtcen = pcc_out[1]*multipoles['Gt0'] 
     Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     plt.rc('font', family='serif', size='13.0')
     ax.text(1,110,'('+str(lab)+')',fontsize = '12')
     ax.plot(r,Gt,'C1--')
     ax.plot(r,Gtcen,'C3')
     ax.plot(r,Gtmiss,'C3--')
     ax.scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
     ax.errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
     ax.set_xscale('log')
     ax.set_yscale('log')
     ax.set_xlabel('r [$h^{-1}$ Mpc]')
     ax.set_ylim(1,200)
     ax.set_xlim(0.3,5)
     ax.axvline(ROUT,ls='--',c='C7')
     ax.xaxis.set_ticks([0.3,1,3])
     ax.set_xticklabels([0.3,1,3])
     ax.yaxis.set_ticks([0.3,10,100])
     ax.set_yticklabels([0.3,10,100])

          


folder    = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/'

f = open(folder+'list_Nbin','r')
lines = f.readlines()

f = open(folder+'list_Mbin','r')
lines = np.append(f.readlines(),lines)

ft, axt = plt.subplots(5,4, figsize=(12,12), sharex=True,sharey=True)
ft.subplots_adjust(hspace=0,wspace=0)

axt[0,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')
axt[1,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')
axt[2,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')
axt[3,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')
axt[4,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')

ax2 = np.reshape(np.array(axt),(20,1))



# lines = lines[:3]+lines[4:]

for j in range(20):
     line = lines[j]
     out = make_plot_misscentred_monopole(line[:-1],folder,ax2[j][0],j+1)
    
plt.savefig(folder+'profiles.pdf',bbox_inches='tight')