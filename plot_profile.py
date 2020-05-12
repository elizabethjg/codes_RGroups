import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits
import corner
from profiles_fit import *

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


def make_plot_misscentred_monopole_tauRmean(file_name,folder):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     file_mcmc = 'monopole_misscentred_'+file_name[:-4]+'out'

     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Rmean   = h['RADIUS_HALO_mean']
     ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))/0.7

     mcmc = (np.loadtxt(folder+file_mcmc)).T
     labels = ['M200','pcc','tau']

     mout    = np.percentile(mcmc[0][1000:], [16, 50, 84])
     pcc_out = np.percentile(mcmc[1][1000:], [16, 50, 84])
     tau_out = np.percentile(mcmc[2][1000:], [16, 50, 84])

     fig = corner.corner(mcmc.T, labels=labels)
     plt.savefig(folder+'plots_monopole_misscentred_tauRmean/corner'+file_mcmc[:-3]+'png')
     
     f, ax = plt.subplots(3, 1, figsize=(6,5))
     ax[0].plot(mcmc[0],'k.',alpha=0.3)
     ax[0].axvline(1000)
     ax[0].axhline(mout[1])
     ax[0].axhline(mout[1] - (mout[1]-mout[0]),ls='--')
     ax[0].axhline(mout[1] + (mout[2]-mout[1]),ls='--')
     ax[1].plot(mcmc[1],'k.',alpha=0.3)
     ax[1].axvline(1000)
     ax[1].axhline(pcc_out[1])
     ax[1].axhline(pcc_out[1] - (pcc_out[1]-pcc_out[0]),ls='--')
     ax[1].axhline(pcc_out[1] + (pcc_out[2]-pcc_out[1]),ls='--')
     ax[2].plot(mcmc[2],'k.',alpha=0.3)
     ax[2].axhline(tau_out[1])
     ax[2].axhline(tau_out[1] - (tau_out[1]-tau_out[0]),ls='--')
     ax[2].axhline(tau_out[1] + (tau_out[2]-tau_out[1]),ls='--')
     ax[2].axvline(1000)
     f.subplots_adjust(hspace=0,wspace=0)
     plt.savefig(folder+'plots_monopole_misscentred_tauRmean/iter'+file_mcmc[:-3]+'png')

     
     
     if h['elM200_NFW'] < 0.:
          eM200 = 10**h['elM200_NFW']
     else:
          eM200 = h['elM200_NFW']
          
     # e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     try:
          print h['N_min'],h['N_max']-1
     except:
          print h['N_GAL_mean']
     
     
     print '####################'
     
     M200   = 10**(mout[1])
     e_M200 = (10**(mout[1])*np.log(10.)*np.diff(mout))
     
     
     
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
     
     
     
     multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                              misscentred = True,s_off = Rmean*tau_out[1],
                              ellip=0,z=zmean,components = ['t'],
                              verbose=False,ncores=2)
     
     Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])     
     
     Gtcen = pcc_out[1]*multipoles['Gt0'] 
     Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     
     f, ax = plt.subplots(2, 1, figsize=(5,8), sharex=True)
     f.subplots_adjust(hspace=0,wspace=0)
     matplotlib.rcParams.update({'font.size': 12})
	
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)

     ax[0].plot(r,Gt,'C1--')
     ax[0].plot(r,Gtcen,'C3')
     ax[0].plot(r,Gtmiss,'C3--')
     ax[0].scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
     ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
     ax[0].set_xscale('log')
     ax[0].set_yscale('log')
     ax[0].set_xlabel('R [Mpc]')
     ax[0].set_ylim(1,200)
     ax[0].set_xlim(0.1,5)
     ax[0].axvline(ROUT,ls='--',c='C7')
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
     plt.savefig(folder+'plots_monopole_misscentred_tauRmean/'+file_name[:-5]+'.png')
     
     
     return Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN, pcc_out[1], np.diff(pcc_out)

def make_plot_misscentred_monopole_tauR200(file_name,folder):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     file_mcmc = 'monopole_3misscentred_'+file_name[:-4]+'out'

     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Rmean   = h['RADIUS_HALO_mean']
     ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))/0.7

     mcmc = (np.loadtxt(folder+file_mcmc)).T
     labels = ['M200','pcc','tau']

     mout    = np.percentile(mcmc[0][1000:], [16, 50, 84])
     pcc_out = np.percentile(mcmc[1][1000:], [16, 50, 84])
     tau_out = np.percentile(mcmc[2][1000:], [16, 50, 84])

     fig = corner.corner(mcmc.T, labels=labels)
     plt.savefig(folder+'plots_monopole_misscentred_tauR200/corner'+file_mcmc[:-3]+'png')

     
     f, ax = plt.subplots(3, 1, figsize=(6,5))
     ax[0].plot(mcmc[0],'k.',alpha=0.3)
     ax[0].axvline(1000)
     ax[0].axhline(mout[1])
     ax[0].axhline(mout[1] - (mout[1]-mout[0]),ls='--')
     ax[0].axhline(mout[1] + (mout[2]-mout[1]),ls='--')
     ax[1].plot(mcmc[1],'k.',alpha=0.3)
     ax[1].axvline(1000)
     ax[1].axhline(pcc_out[1])
     ax[1].axhline(pcc_out[1] - (pcc_out[1]-pcc_out[0]),ls='--')
     ax[1].axhline(pcc_out[1] + (pcc_out[2]-pcc_out[1]),ls='--')
     ax[2].plot(mcmc[2],'k.',alpha=0.3)
     ax[2].axhline(tau_out[1])
     ax[2].axhline(tau_out[1] - (tau_out[1]-tau_out[0]),ls='--')
     ax[2].axhline(tau_out[1] + (tau_out[2]-tau_out[1]),ls='--')
     ax[2].axvline(1000)
     f.subplots_adjust(hspace=0,wspace=0)
     plt.savefig(folder+'plots_monopole_misscentred_tauR200/iter'+file_mcmc[:-3]+'png')

     
     
     if h['elM200_NFW'] < 0.:
          eM200 = 10**h['elM200_NFW']
     else:
          eM200 = h['elM200_NFW']
          
     # e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     try:
          print h['N_min'],h['N_max']-1
     except:
          print h['N_GAL_mean']
     
     
     print '####################'
     
     M200   = 10**(mout[1])
     e_M200 = (10**(mout[1])*np.log(10.)*np.diff(mout))
     
     # Compute cosmological parameters
     cosmo = LambdaCDM(H0=0.7*100, Om0=0.3, Ode0=0.7)
     H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
     roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
     roc_mpc  = roc*((pc*1.0e6)**3.0)

     R200 = r200_nfw(M200,roc_mpc)
     
     
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
     
     
     
     multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                              misscentred = True,s_off = R200*tau_out[1],
                              ellip=0,z=zmean,components = ['t'],
                              verbose=False,ncores=2)
     
     Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])     
     
     Gtcen = pcc_out[1]*multipoles['Gt0'] 
     Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     
     f, ax = plt.subplots(2, 1, figsize=(5,8), sharex=True)
     f.subplots_adjust(hspace=0,wspace=0)
     matplotlib.rcParams.update({'font.size': 12})
	
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)

     ax[0].plot(r,Gt,'C1--')
     ax[0].plot(r,Gtcen,'C3')
     ax[0].plot(r,Gtmiss,'C3--')
     ax[0].scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
     ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
     ax[0].set_xscale('log')
     ax[0].set_yscale('log')
     ax[0].set_xlabel('R [Mpc]')
     ax[0].set_ylim(1,200)
     ax[0].set_xlim(0.1,5)
     ax[0].axvline(ROUT,ls='--',c='C7')
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
     plt.savefig(folder+'plots_monopole_misscentred_tauR200/'+file_name[:-5]+'.png')
     
     return Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN, pcc_out[1], np.diff(pcc_out)

def make_plot_misscentred_monopole_soff(file_name,folder):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     file_mcmc = 'monopole_4misscentred_'+file_name[:-4]+'out'

     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Rmean   = h['RADIUS_HALO_mean']
     ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))/0.7

     mcmc = (np.loadtxt(folder+file_mcmc)).T
     labels = ['M200','pcc','s_off']

     mout    = np.percentile(mcmc[0][1000:], [16, 50, 84])
     pcc_out = np.percentile(mcmc[1][1000:], [16, 50, 84])
     soff_out = np.percentile(mcmc[2][1000:], [16, 50, 84])

     fig = corner.corner(mcmc.T, labels=labels)
     plt.savefig(folder+'plots_monopole_misscentred_soff/corner'+file_mcmc[:-3]+'png')
   
     f, ax = plt.subplots(3, 1, figsize=(6,5))
     ax[0].plot(mcmc[0],'k.',alpha=0.3)
     ax[0].axvline(1000)
     ax[0].axhline(mout[1])
     ax[0].axhline(mout[1] - (mout[1]-mout[0]),ls='--')
     ax[0].axhline(mout[1] + (mout[2]-mout[1]),ls='--')
     ax[1].plot(mcmc[1],'k.',alpha=0.3)
     ax[1].axvline(1000)
     ax[1].axhline(pcc_out[1])
     ax[1].axhline(pcc_out[1] - (pcc_out[1]-pcc_out[0]),ls='--')
     ax[1].axhline(pcc_out[1] + (pcc_out[2]-pcc_out[1]),ls='--')
     ax[2].plot(mcmc[2],'k.',alpha=0.3)
     ax[2].axhline(soff_out[1])
     ax[2].axhline(soff_out[1] - (soff_out[1]-soff_out[0]),ls='--')
     ax[2].axhline(soff_out[1] + (soff_out[2]-soff_out[1]),ls='--')
     ax[2].axvline(1000)
     f.subplots_adjust(hspace=0,wspace=0)
     plt.savefig(folder+'plots_monopole_misscentred_soff/iter'+file_mcmc[:-3]+'png')

     
     
     if h['elM200_NFW'] < 0.:
          eM200 = 10**h['elM200_NFW']
     else:
          eM200 = h['elM200_NFW']
          
     # e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     try:
          print h['N_min'],h['N_max']-1
     except:
          print h['N_GAL_mean']
     
     
     print '####################'
     
     M200   = 10**(mout[1])
     e_M200 = (10**(mout[1])*np.log(10.)*np.diff(mout))
     
     
     
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
     
     
     
     multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                              misscentred = True,s_off = soff_out[1],
                              ellip=0,z=zmean,components = ['t'],
                              verbose=False,ncores=2)
     
     Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])     
     
     Gtcen = pcc_out[1]*multipoles['Gt0'] 
     Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     
     f, ax = plt.subplots(2, 1, figsize=(5,8), sharex=True)
     f.subplots_adjust(hspace=0,wspace=0)
     matplotlib.rcParams.update({'font.size': 12})
	
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)

     ax[0].plot(r,Gt,'C1--')
     ax[0].plot(r,Gtcen,'C3')
     ax[0].plot(r,Gtmiss,'C3--')
     ax[0].scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
     ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
     ax[0].set_xscale('log')
     ax[0].set_yscale('log')
     ax[0].set_xlabel('R [Mpc]')
     ax[0].set_ylim(1,200)
     ax[0].set_xlim(0.1,5)
     ax[0].axvline(ROUT,ls='--',c='C7')
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
     plt.savefig(folder+'plots_monopole_misscentred_soff/'+file_name[:-5]+'.png')
     
     return Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN, pcc_out[1], np.diff(pcc_out)


def make_plot_misscentred_monopole_pcc(file_name,folder,plot = False):
          

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
     
     mcmc = (np.loadtxt(folder+file_mcmc)).T
     labels = ['M200','pcc']

     mout    = np.percentile(mcmc[0][1000:], [16, 50, 84])
     pcc_out = np.percentile(mcmc[1][1000:], [16, 50, 84])
     
     
     if h['elM200_NFW'] < 0.:
          eM200 = 10**h['elM200_NFW']
     else:
          eM200 = h['elM200_NFW']
          
     try:
          MDYN = 10**h['lMASS_DYN_mean']
     except:
          MDYN = -1.
          
     # e_M200 = (M200*eM200)/np.log(10.)
          
     zmean   = h['Z_MEAN']
     
     print '####################'
     print file_name
     print 'lM200 ',h['lM200_NFW']
     
     try:
          print h['N_min'],h['N_max']-1
          nmin, nmax = h['N_min'],h['N_max']-1
     except:
          print h['N_GAL_mean']
          nmin = h['N_GAL_mean']
          nmax = h['N_GAL_mean']
     
     print 'C_BG, min, mas'
     print h['C_BG_min'],h['C_BG_max']
     print 'lMH, min, mas'
     print h['lMH_min'],h['lMH_max']
     
     
     print '####################'
     
     M200   = 10**(mout[1])
     e_M200 = (10**(mout[1])*np.log(10.)*np.diff(mout))
     
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)
     
     if plot:
     
          r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
          
          multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                                   misscentred = True,s_off = soff,
                                   ellip=0,z=zmean,components = ['t'],
                                   verbose=False,ncores=2)
          
          Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])     
          
          Gtcen = pcc_out[1]*multipoles['Gt0'] 
          Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     
          fig = corner.corner(mcmc.T, labels=labels)
          plt.savefig(folder+'plots_monopole_misscentred_pcconly/corner'+file_mcmc[:-3]+'png')
          
          f, ax = plt.subplots(2, 1, figsize=(6,3))
          ax[0].plot(mcmc[0],'k.',alpha=0.3)
          ax[0].axvline(1000)
          ax[0].axhline(mout[1])
          ax[0].axhline(mout[1] - (mout[1]-mout[0]),ls='--')
          ax[0].axhline(mout[1] + (mout[2]-mout[1]),ls='--')
          ax[1].plot(mcmc[1],'k.',alpha=0.3)
          ax[1].axvline(1000)
          ax[1].axhline(pcc_out[1])
          ax[1].axhline(pcc_out[1] - (pcc_out[1]-pcc_out[0]),ls='--')
          ax[1].axhline(pcc_out[1] + (pcc_out[2]-pcc_out[1]),ls='--')
          f.subplots_adjust(hspace=0,wspace=0)
          plt.savefig(folder+'plots_monopole_misscentred_pcconly/iter'+file_mcmc[:-3]+'png')
          
          f, ax = plt.subplots(2, 1, figsize=(5,8), sharex=True)
          f.subplots_adjust(hspace=0,wspace=0)
          matplotlib.rcParams.update({'font.size': 12})
     
          ax[0].plot(r,Gt,'C1--')
          ax[0].plot(r,Gtcen,'C3')
          ax[0].plot(r,Gtmiss,'C3--')
          ax[0].scatter(p.Rp,p.DSigma_T,facecolor='none',edgecolors='0.4')
          ax[0].errorbar(p.Rp,p.DSigma_T,yerr=p.error_DSigma_T,fmt = 'none',ecolor='0.4')
          ax[0].set_xscale('log')
          ax[0].set_yscale('log')
          ax[0].set_xlabel('R [Mpc]')
          ax[0].set_ylim(1,200)
          ax[0].set_xlim(0.1,5)
          ax[0].axvline(ROUT,ls='--',c='C7')
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
          plt.savefig(folder+'plots_monopole_misscentred_pcconly/'+file_name[:-5]+'.png')
     
     return Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN, pcc_out[1], np.diff(pcc_out), MDYN/1.e14, nmin, nmax, zmean


folder = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_FOF/'
     
f = open(folder+'list_Mbin','r')
# f = open(folder+'list_m1','r')
lines = f.readlines()

MH    = np.array([])
MDYN    = np.array([])
MNFW  = np.array([])
eMNFW = np.array([])
N_mean = np.array([])
Nlenses = np.array([])
SN = np.array([])
pcc = np.array([])
e_pcc = np.array([])
NMIN = np.array([])
NMAX = np.array([])
ZMEAN = np.array([])
# lines = lines[:3]+lines[4:]

for line in lines:
     
     Mhalo, M200, eM200, Nmean, Nlens, sn, pcc_, e_pcc_,mdyn,nmin,nmax,zmean = make_plot_misscentred_monopole_pcc(line[:-1],folder,True)
     
     MH    = np.append(MH,Mhalo)
     MNFW  = np.append(MNFW,M200*0.7)
     eMNFW = np.append(eMNFW,[eM200*0.7])
     N_mean = np.append(N_mean,Nmean)
     Nlenses = np.append(Nlenses,Nlens)
     SN      = np.append(SN,sn)
     pcc     = np.append(pcc,pcc_)
     e_pcc   = np.append(e_pcc,e_pcc_)
     MDYN   = np.append(MDYN,mdyn)
     NMIN   = np.append(NMIN,nmin)
     NMAX   = np.append(NMAX,nmax)
     ZMEAN   = np.append(ZMEAN,zmean)

eMNFW = np.reshape(eMNFW,(len(lines),2))
e_pcc = np.reshape(e_pcc,(len(lines),2))

lMNFW = np.log10(MNFW*1.e14)
lMDYN = np.log10(MDYN*1.e14)
lMH   = np.log10(MH*1.e14)
elMNFW = [eMNFW.T[0]/(MNFW*np.log(10.)),eMNFW.T[1]/(MNFW*np.log(10.))]

j = np.argsort(N_mean)

out = np.array([NMIN[j],NMAX[j],N_mean[j],ZMEAN[j],lMH[j],lMDYN[j],
                lMNFW[j],elMNFW[0][j],elMNFW[1][j],pcc[j],
                e_pcc.T[0][j],e_pcc.T[1][j],Nlenses[j]])


f1=open(folder+'Lens_results_Mbin.out','w')
f1.write('# NMIN     NMAX    Nmean  zmean  logMHALO    logMDYN   logMNFW   elogMNFW-   elogMNFW+   pcc   epcc-   epcc+   Nlenses \n')
np.savetxt(f1,out.T,fmt = ['%4i']*2+['%12.4f']*10+['%6i'])
f1.close()
