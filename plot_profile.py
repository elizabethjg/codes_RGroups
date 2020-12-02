import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/home/elizabeth/Documentos/posdoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/eli/python_codes')
sys.path.append('/home/elizabeth/python_codes')
sys.path.append('/mnt/clemente/lensing/python_codes')
import numpy as np
from matplotlib import *
from multipoles_shear import *
from astropy.io import fits
from profiles_fit import *
import os
import argparse
import corner

def make_plot_centred_monopole(file_name,folder,samples,ncores):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     M200    = 10**h['lM200_NFW']
     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Vdisp   = h['VDISP_DYN_mean']
     
     eM200 = 10**h['elM200_NFW']          
     try:
          MDYN = 10**h['lMASS_DYN_mean']
     except:
          MDYN = -1.
          
     e_M200 = (M200*eM200)/np.log(10.)
          
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
     
     cmin, cmax = h['C_BG_min'],h['C_BG_max']
     Mmin, Mmax = h['lMH_min'],h['lMH_max']
     zmin, zmax = h['z_min'],h['z_max']
     
     print 'C_BG, min, max'
     print cmin, cmax
     print 'lMH, min, max'
     print Mmin, Mmax
     print 'z, min, max'
     print zmin, zmax

     
     print '####################'

     
     r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
     
     
     
     out = multipole_shear_parallel(r,M200=M200,z=zmean,
                                   ellip=0.,misscentred=False,
                                   ncores=ncores)
     
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
     ax[0].set_xlim(0.3,5)
     ax[0].xaxis.set_ticks([0.3,1,5])
     ax[0].set_xticklabels([0.3,1,5])
     ax[0].yaxis.set_ticks([0.3,10,100])
     ax[0].set_yticklabels([0.3,10,100])
     ax[0].set_ylabel(r'$\Delta \Sigma_T [h_{70}M_\odot\,\rm{pc}^{-2}]$')
     
     ax[1].plot([0,5],[0,0],'k--')
     ax[1].scatter(p.Rp,p.DSigma_X,facecolor='none',edgecolors='0.4')
     ax[1].errorbar(p.Rp,p.DSigma_X,yerr=p.error_DSigma_X,fmt = 'none',ecolor='0.4')
     ax[1].set_xscale('log')
     ax[1].set_xlabel('R [mpc]')
     ax[1].set_ylim(-50,50)
     ax[1].set_xlim(0.3,5)
     ax[1].xaxis.set_ticks([0.3,1,5])
     ax[1].set_xticklabels([0.3,1,5])
     ax[1].yaxis.set_ticks([-25,0,25])
     ax[1].set_yticklabels([-25,0,25])
     ax[1].set_ylabel(r'$\Delta \Sigma_\times [h_{70}M_\odot\,\rm{pc}^{-2}]$')
     matplotlib.rcParams.update({'font.size': 12})
     plt.savefig(folder+'plots_'+samples+'/'+file_name[:-5]+'_cen.png')
     
     return [Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN, 
            1., np.array([0,0]), MDYN/1.e14, nmin, nmax, 
            Mmin, Mmax, zmin, zmax, zmean, 1.0, Vdisp]



def make_plot_misscentred_monopole(file_name,folder,samples,plot = False,ymiss = False,ncores = 32):
          

     profile = fits.open(folder+file_name)
     h       = profile[1].header
     p       = profile[1].data

     if ymiss:
          file_mcmc = 'monopole_pcconly_ymiss_'+file_name[:-4]+'out'
          tau       = 0.3
     else:
          file_mcmc = 'monopole_pcconly_boost_'+file_name[:-4]+'out'
          tau       = 0.4
          
     Mhalo   = 10**h['lMASS_HALO_mean']
     Nmean   = h['N_GAL_mean']
     Nlens   = h['N_LENSES']
     Rmean   = h['RADIUS_HALO_mean']
     ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))
     soff = tau*Rmean
     Vdisp   = h['VDISP_DYN_mean']
     
     try:
          mcmc = (np.loadtxt(folder+file_mcmc)).T
     except:
          return 'out' 
          
     labels = ['M200','pcc']

     mout    = np.percentile(mcmc[0][1000:], [16, 50, 84])
     pcc_out = np.percentile(mcmc[1][1000:], [16, 50, 84])
     
     
     eM200 = 10**h['elM200_NFW']          
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
     
     SN = np.mean(p.DSigma_T/p.error_DSigma_T)

	######################
	#COMPUTE CHI_SQUARE

     mr = p.Rp < ROUT
     multipoles = multipole_shear_parallel(p.Rp[mr],M200=10**mout[1],
                                   misscentred = True,s_off = soff,h = 1.0,
                                   ellip=0,z=zmean,components = ['t'],
                                   verbose=False,ncores=ncores,Yanmiss=ymiss)
     
     modelt_t = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc_out[1])
     
     chi_t = chi_red(modelt_t,p.DSigma_T[mr],p.error_DSigma_T[mr],2)

     print chi_t
     print '####################'
	#--------------------
    
     if plot:

          print 'PLOTTING...'
          
          if ymiss:
               folder_plot = 'plots_ymiss_'+samples+'/'
          else:
               folder_plot = 'plots_'+samples+'/'

          print folder+folder_plot
          print '------------'

          
          r  = np.logspace(np.log10(0.05),np.log10(5.5),20)
          
          multipoles = multipole_shear_parallel(r,M200=10**mout[1],
                                   misscentred = True,s_off = soff,h = 1.0,
                                   ellip=0,z=zmean,components = ['t'],
                                   verbose=False,ncores=ncores,Yanmiss=ymiss)
          
          Gt    = model_Gamma(multipoles,'t',misscentred=True,pcc=pcc_out[1])     
          
          Gtcen = pcc_out[1]*multipoles['Gt0'] 
          Gtmiss = (1-pcc_out[1])*multipoles['Gt_off']
     
          fig = corner.corner(mcmc.T, labels=labels)
          plt.savefig(folder+folder_plot+'corner'+file_mcmc[:-3]+'png')
          
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
          plt.savefig(folder+folder_plot+file_mcmc[:-3]+'png')
          
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
          ax[0].set_xlim(0.3,5)
          ax[0].axvline(ROUT,ls='--',c='C7')
          ax[0].xaxis.set_ticks([0.3,1,5])
          ax[0].set_xticklabels([0.3,1,5])
          ax[0].yaxis.set_ticks([0.3,10,100])
          ax[0].set_yticklabels([0.3,10,100])
          ax[0].set_ylabel(r'$\Delta \Sigma_T [h_{70}M_\odot\,\rm{pc}^{-2}]$')
          
          ax[1].plot([0,5],[0,0],'k--')
          ax[1].scatter(p.Rp,p.DSigma_X,facecolor='none',edgecolors='0.4')
          ax[1].errorbar(p.Rp,p.DSigma_X,yerr=p.error_DSigma_X,fmt = 'none',ecolor='0.4')
          ax[1].set_xscale('log')
          ax[1].set_xlabel('R [mpc]')
          ax[1].set_ylim(-50,50)
          ax[1].set_xlim(0.3,5)
          ax[1].xaxis.set_ticks([0.3,1,5])
          ax[1].set_xticklabels([0.3,1,5])
          ax[1].yaxis.set_ticks([-25,0,25])
          ax[1].set_yticklabels([-25,0,25])
          ax[1].set_ylabel(r'$\Delta \Sigma_\times [h_{70}M_\odot\,\rm{pc}^{-2}]$')
          matplotlib.rcParams.update({'font.size': 12})
          plt.savefig(folder+folder_plot+file_name[:-5]+'_miss.png')
     
     return [Mhalo/1.e14, M200/1.e14, e_M200/1.e14, Nmean, Nlens, SN, 
            pcc_out[1], np.diff(pcc_out), MDYN/1.e14, nmin, nmax, 
            Mmin, Mmax, zmin, zmax, zmean, chi_t, Vdisp]


parser = argparse.ArgumentParser()
parser.add_argument('-sample', action='store', dest='sample',default='Mbin')
parser.add_argument('-ymiss', action='store', dest='ymiss', default='False')
parser.add_argument('-plot', action='store', dest='plot', default='True')
parser.add_argument('-only_central', action='store', dest='cen', default='False')
parser.add_argument('-ncores', action='store', dest='ncores', default=3)
args = parser.parse_args()

# folder    = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/'
folder    = '/home/elizabeth/Documentos/posdoc/Rgroups/profiles_indcat/'
# folder    = '/mnt/clemente/lensing/RodriguezGroups/N_all/'
samples   = args.sample
ncores    = args.ncores
if 'True' in args.ymiss:
     ymiss = True
elif 'False' in args.ymiss:
     ymiss = False

if 'True' in args.plot:
     makeplots = True
elif 'False' in args.plot:
     makeplots = False

if 'True' in args.cen:
     cen = True
elif 'False' in args.cen:
     cen = False


print folder
print 'sample', samples
print 'ncores',ncores
print 'ymiss', ymiss
print 'plot', makeplots


f = open(folder+'list_'+samples,'r')
# f = open(folder+'list_m1','r')
lines = f.readlines()

if ymiss:
     print 'Y miss!!!!!!!!!!!!!!!!!'
     os.system('mkdir '+folder+'plots_ymiss_'+samples)
else:
     os.system('mkdir '+folder+'plots_'+samples)


MNFW  = np.array([])
eMNFW = np.array([])

Nlenses = np.array([])
SN = np.array([])
pcc = np.array([])
e_pcc = np.array([])

NMIN = np.array([])
NMAX = np.array([])
MMIN = np.array([])
MMAX = np.array([])
ZMIN = np.array([])
ZMAX = np.array([])

Zmean = np.array([])
Mmean = np.array([])
Nmean = np.array([])
MDYN  = np.array([])
Vdisp = np.array([])

CHI2  = np.array([])

# lines = lines[:3]+lines[4:]

nbins = 0

for line in lines:
     
     
     try:
          if cen:
               out = make_plot_centred_monopole(line[:-1],folder,samples,ncores)
          else:
               out = make_plot_misscentred_monopole(line[:-1],folder,samples,makeplots,ymiss,ncores)
     except:
          continue
    
     nbins += 1
     Mmean    = np.append(Mmean,out[0])
     MNFW     = np.append(MNFW,out[1])
     eMNFW = np.append(eMNFW,[out[2]])
     
     Nmean = np.append(Nmean,out[3])
     Nlenses = np.append(Nlenses,out[4])
     SN      = np.append(SN,out[5])
     pcc     = np.append(pcc,out[6])
     e_pcc   = np.append(e_pcc,out[7])
     MDYN   = np.append(MDYN,out[8])
     
     NMIN   = np.append(NMIN,out[9])
     NMAX   = np.append(NMAX,out[10])
     MMIN   = np.append(MMIN,10**out[11]/1.e14)
     MMAX   = np.append(MMAX,10**out[12]/1.e14)
     ZMIN   = np.append(ZMIN,out[13])
     ZMAX   = np.append(ZMAX,out[14])
     
     Zmean   = np.append(Zmean,out[15])
     CHI2    = np.append(CHI2,out[16])
     Vdisp   = np.append(Vdisp,out[17])


eMNFW = np.reshape(eMNFW,(nbins,2)).T
e_pcc = np.reshape(e_pcc,(nbins,2)).T

#lMNFW = np.log10(MNFW*1.e14)
#lMDYN = np.log10(MDYN*1.e14)
#lMH   = np.log10(MH*1.e14)
#elMNFW = [eMNFW.T[0]/(MNFW*np.log(10.)),eMNFW.T[1]/(MNFW*np.log(10.))]

j = np.argsort(Nmean)

out = np.array([Nlenses[j], NMIN[j],NMAX[j],Nmean[j],
                MMIN[j],MMAX[j],Mmean[j],
                ZMIN[j],ZMAX[j],Zmean[j],
                MDYN[j], MNFW[j], eMNFW[0][j],eMNFW[1][j],
                pcc[j], e_pcc[0][j],e_pcc[1][j],CHI2[j],Vdisp[j]])

if ymiss:
     f1=open(folder+'Lens_'+samples+'_ymiss.out','w')
else:
     f1=open(folder+'Lens_'+samples+'.out','w')
     
f1.write('# Nlenses    Nmin    Nmax   Nmean   Mmin        Mmax         Mmean       zmin         zmax         zmean         Mdyn         Mlens       eMlens                    pcc          e_pcc                        chi2       Vdisp\n')
np.savetxt(f1,out.T,fmt = ['%4i']*3+['%12.4f']*16)
f1.close()
