import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits

def make_plot_centred_monopole(file_name,folder):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     M200    = 10**h['lM200_NFW']
     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
    
     if h['elM200_NFW'] < 0.:
          eM200 = 10**h['elM200_NFW']
     else:
          eM200 = h['elM200_NFW']
          
     e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     try:
          print h['N_min'],h['N_max']-1
     except:
          print h['N_GAL_mean']

     
     print '####################'

     
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
     
     
     
     out = multipole_shear_parallel(r,M200=M200,z=zmean,
                                   ellip=0.,misscentred=False,
                                   ncores=2)
     
     Gt    = model_Gamma(out,'t',misscentred=False)     
     
     f, ax = plt.subplots(2, 1, figsize=(5,8), sharex=True)
     f.subplots_adjust(hspace=0,wspace=0)
     matplotlib.rcParams.update({'font.size': 12})
	
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)

     ax[0].plot(r,Gt,'C1--')
     ax[0].scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
     ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
     ax[0].set_xscale('log')
     ax[0].set_yscale('log')
     ax[0].set_xlabel('R [Mpc]')
     ax[0].set_ylim(1,200)
     ax[0].set_xlim(0.1,5)
     ax[0].xaxis.set_ticks([0.1,1,5])
     ax[0].set_xticklabels([0.1,1,5])
     ax[0].yaxis.set_ticks([0.1,10,100])
     ax[0].set_yticklabels([0.1,10,100])
     ax[0].set_ylabel(r'$\Delta \Sigma_T [h_{70}M_\odot\,\rm{pc}^{-2}]$')
     
     ax[1].plot([0,5],[0,0],'k--')
     ax[1].scatter(p.Rp,p.DSigma_X,facecolor='none',edgecolors='0.4')
     ax[1].errorbar(p.Rp,p.DSigma_X,yerr=p.error_DSigma_X,fmt = 'none',ecolor='0.4')
     ax[1].set_xscale('log')
     ax[1].set_xlabel('R [mpc]')
     ax[1].set_ylim(-50,50)
     ax[1].set_xlim(0.1,5)
     ax[1].xaxis.set_ticks([0.1,1,5])
     ax[1].set_xticklabels([0.1,1,5])
     ax[1].yaxis.set_ticks([-25,0,25])
     ax[1].set_yticklabels([-25,0,25])
     ax[1].set_ylabel(r'$\Delta \Sigma_\times [h_{70}M_\odot\,\rm{pc}^{-2}]$')
     matplotlib.rcParams.update({'font.size': 12})
     plt.savefig(folder+'plots_monopole_centred/'+file_name[:-5]+'.png')
     
     return Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN

def make_plot_misscentred_monopole(file_name,folder):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     M200    = 10**h['lM200_NFW']
     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
    
     if h['elM200_NFW'] < 0.:
          eM200 = 10**h['elM200_NFW']
     else:
          eM200 = h['elM200_NFW']
          
     e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     try:
          print h['N_min'],h['N_max']-1
     except:
          print h['N_GAL_mean']

     
     print '####################'

     
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
     
     
     
     out = multipole_shear_parallel(r,M200=M200,z=zmean,
                                   ellip=0.,misscentred=False,
                                   ncores=2)
     
     Gt    = model_Gamma(out,'t',misscentred=False)     
     
     f, ax = plt.subplots(2, 1, figsize=(5,8), sharex=True)
     f.subplots_adjust(hspace=0,wspace=0)
     matplotlib.rcParams.update({'font.size': 12})
	
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)

     ax[0].plot(r,Gt,'C1--')
     ax[0].scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
     ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
     ax[0].set_xscale('log')
     ax[0].set_yscale('log')
     ax[0].set_xlabel('R [Mpc]')
     ax[0].set_ylim(1,200)
     ax[0].set_xlim(0.1,5)
     ax[0].xaxis.set_ticks([0.1,1,5])
     ax[0].set_xticklabels([0.1,1,5])
     ax[0].yaxis.set_ticks([0.1,10,100])
     ax[0].set_yticklabels([0.1,10,100])
     ax[0].set_ylabel(r'$\Delta \Sigma_T [h_{70}M_\odot\,\rm{pc}^{-2}]$')
     
     ax[1].plot([0,5],[0,0],'k--')
     ax[1].scatter(p.Rp,p.DSigma_X,facecolor='none',edgecolors='0.4')
     ax[1].errorbar(p.Rp,p.DSigma_X,yerr=p.error_DSigma_X,fmt = 'none',ecolor='0.4')
     ax[1].set_xscale('log')
     ax[1].set_xlabel('R [mpc]')
     ax[1].set_ylim(-50,50)
     ax[1].set_xlim(0.1,5)
     ax[1].xaxis.set_ticks([0.1,1,5])
     ax[1].set_xticklabels([0.1,1,5])
     ax[1].yaxis.set_ticks([-25,0,25])
     ax[1].set_yticklabels([-25,0,25])
     ax[1].set_ylabel(r'$\Delta \Sigma_\times [h_{70}M_\odot\,\rm{pc}^{-2}]$')
     matplotlib.rcParams.update({'font.size': 12})
     plt.savefig(folder+'plots_monopole_centred/'+file_name[:-5]+'.png')
     
     return Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN



folder = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/'
     
f = open(folder+'list_names','r')
lines = f.readlines()

MH    = np.array([])
MNFW  = np.array([])
eMNFW = np.array([])
N_mean = np.array([])
Nlenses = np.array([])
SN = np.array([])

for line in lines:
     
     Mhalo, M200, eM200, Nmean, Nlens, sn = make_plot_centred_monopole(line[:-1],folder)
     
     MH    = np.append(MH,Mhalo)
     MNFW  = np.append(MNFW,M200*0.7)
     eMNFW = np.append(eMNFW,eM200*0.7)
     N_mean = np.append(N_mean,Nmean)
     Nlenses = np.append(Nlenses,Nlens)
     SN      = np.append(SN,sn)


lMNFW = np.log10(MNFW*1.e14)
lMH   = np.log10(MH*1.e14)
elMNFW = eMNFW/(MNFW*np.log(10.))
