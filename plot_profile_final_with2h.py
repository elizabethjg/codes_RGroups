import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/Documentos/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/python_codes')
import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits
import corner
from profiles_fit import *
import os
import argparse



def make_plot_misscentred_monopole(file_name,folder,at,lab):
          
     
     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data
     
     file_mcmc = 'monopole_pcconly_'+file_name[:-4]+'out'
     file_2h   = '../second_halo_RodriguezFOF/smooth_second_halo_'+str(lab)+'.dat'
          
     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Rmean   = h['RADIUS_HALO_mean']
     ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))
     soff = 0.4*Rmean
     Vdisp   = h['VDISP_DYN_mean']
     
     mcmc = (np.loadtxt(folder+file_mcmc)).T
     two_halo = (np.loadtxt(file_2h)).T

     bines = np.logspace(np.log10(300.),np.log10(5000.),num=len(p)+1)
     
     area = np.pi*np.diff(bines**2)
     ngal = p.NGAL_w
     
     d = ngal/area
     
     fcl = ((d - d[-1])*area)/ngal
     
     bcorr = 1./(1.-fcl)
     
     bcorr[bcorr < 1.] = 1.
     
     print 'Correcting for boost'
     print bcorr
          
     
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
          
     r  = two_halo[0]
          
     multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                                   misscentred = True,s_off = soff,h=1,
                                   ellip=0,z=zmean,components = ['t'],
                                   verbose=False,ncores=30)
          
     Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])
     Gt2h  = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1]) + two_halo[1]   
          
     Gtcen = pcc_out[1]*multipoles['Gt0'] 
     Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     plt.rc('font', family='serif', size='13.0')
     at.text(1,110,'['+str('%.1f' % Mmin)+','+str('%.1f' % Mmax)+')',fontsize = '12')
     at.plot(r,Gt2h,'C1')
     at.plot(r,Gt,'C1--')
     at.plot(r,Gtcen,'C3')
     at.plot(r,Gtmiss,'C3--')
     at.plot(r,two_halo[1],'C3:')
     if lab < 8:
          at.scatter(p.Rp,p.DSigma_T*bcorr,facecolor='none',edgecolors='C7')
     else:
          at.plot(p.Rp,p.DSigma_T,'C7o')
     at.errorbar(p.Rp,p.DSigma_T*bcorr,yerr=p.error_DSigma_T*bcorr,fmt = 'none',ecolor='0.4')
     at.plot([ROUT,ROUT],[0.1,80],'C7--')
     at.set_xscale('log')
     at.set_yscale('log')
     at.set_xlabel('r [$h^{-1}$ Mpc]')
     at.set_ylim(0.1,200)
     at.set_xlim(0.3,5)
     at.xaxis.set_ticks([0.4,1,3])
     at.set_xticklabels([0.4,1,3])
     at.yaxis.set_ticks([0.1,1,10,100])
     at.set_yticklabels([0.1,1,10,100])


          
# folder    = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/'
folder    = '../profiles_indcat/'
# folder    = '/mnt/clemente/lensing/RodriguezGroups/N_all/'

# f = open(folder+'list_Mbinb','r')
# lines = f.readlines()
sample = 'Mbin'
f = open(folder+'list_'+sample,'r')
lines = f.readlines()

ft, axt = plt.subplots(2,4, figsize=(12,6), sharey=True)
ft.subplots_adjust(hspace=0,wspace=0)

axt[0,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')
axt[1,0].set_ylabel('$\Delta \Sigma$ [$h M_\odot$ /pc$^{2}$]')

at2 = np.reshape(np.array(axt),(8,1))

# lines = lines[:3]+lines[4:]

for j in range(len(lines)):
     line = lines[j]
     out = make_plot_misscentred_monopole(line[:-1],folder,at2[j][0],j)

axt[1,3].axis('off')
axt[1,3].xaxis.set_ticks([0.4,1,3])
axt[1,3].set_xticklabels([0.4,1,3])


    
ft.savefig(folder+'profiles_'+sample+'_2h.pdf',bbox_inches='tight')

