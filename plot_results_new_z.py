import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from matplotlib import rc
from astropy.cosmology import LambdaCDM
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 16})
plt.rc('font', family='serif', size='16.0')

cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)
        
#parameters
cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)


def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er


pcc_s   = (np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/results_simu/pcc.txt').T)
masa_s   = (np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/results_simu/puntos-masamockmasaab_mass.txt').T)
z_s   = (np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/results_simu/puntos-masamockmasaab_z.txt').T)

lMH_simu = pcc_s[2]
pcc_simu = 1-(pcc_s[-1]/100.)

out_T   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_Mbin.out').T
out_H    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zH.out').T
out_L    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zL.out').T

MH_H   = out_H[6]
MH_L   = out_L[6]
MH_T = out_T[6]

lMH_H  = np.log10(out_H[6]*1.e14)
lMH_L  = np.log10(out_L[6]*1.e14)
lMH_T  = np.log10(out_T[6]*1.e14)

M200_H  = out_H[11] 
M200_L  = out_L[11] 
M200_T  = out_T[11] 

eM200_T   = np.array([out_T[12],out_T[13]])
eM200_H   = np.array([out_H[12],out_H[13]])
eM200_L   = np.array([out_L[12],out_L[13]])

lM200_H  = np.log10(out_H[11]*1.e14)
lM200_L  = np.log10(out_L[11]*1.e14)
lM200_T  = np.log10(out_T[11]*1.e14)

elM200_T = np.array([out_T[12]/(out_T[11]*np.log(10.)),out_T[13]/(out_T[11]*np.log(10.))])
elM200_H = np.array([out_H[12]/(out_H[11]*np.log(10.)),out_H[13]/(out_H[11]*np.log(10.))])
elM200_L = np.array([out_L[12]/(out_L[11]*np.log(10.)),out_L[13]/(out_L[11]*np.log(10.))])

pcc_T = out_T[14]
pcc_H = out_H[14]
pcc_L = out_L[14]

epcc_T     = np.array([out_T[15],out_T[16]])
epcc_H     = np.array([out_H[15],out_H[16]])
epcc_L     = np.array([out_L[15],out_L[16]])


########## CON SAMPLES ###########


out_TcM   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_Mbin_cM.out').T
out_HcM    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zH_cM.out').T
out_LcM    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zL_cM.out').T

MH_HcM   = out_HcM[6]
MH_LcM   = out_LcM[6]
MH_TcM = out_TcM[6]

lMH_HcM  = np.log10(out_HcM[6]*1.e14)
lMH_LcM  = np.log10(out_LcM[6]*1.e14)
lMH_TcM  = np.log10(out_TcM[6]*1.e14)

M200_HcM  = out_HcM[11] 
M200_LcM  = out_LcM[11] 
M200_TcM  = out_TcM[11] 

eM200_TcM   = np.array([out_TcM[12],out_TcM[13]])
eM200_HcM   = np.array([out_HcM[12],out_HcM[13]])
eM200_LcM   = np.array([out_LcM[12],out_LcM[13]])

lM200_HcM  = np.log10(out_HcM[11]*1.e14)
lM200_LcM  = np.log10(out_LcM[11]*1.e14)
lM200_TcM  = np.log10(out_TcM[11]*1.e14)

elM200_TcM = np.array([out_TcM[12]/(out_TcM[11]*np.log(10.)),out_TcM[13]/(out_TcM[11]*np.log(10.))])
elM200_HcM = np.array([out_HcM[12]/(out_HcM[11]*np.log(10.)),out_HcM[13]/(out_HcM[11]*np.log(10.))])
elM200_LcM = np.array([out_LcM[12]/(out_LcM[11]*np.log(10.)),out_LcM[13]/(out_LcM[11]*np.log(10.))])

pcc_TcM = out_TcM[14]
pcc_HcM = out_HcM[14]
pcc_LcM = out_LcM[14]

epcc_TcM     = np.array([out_TcM[15],out_TcM[16]])
epcc_HcM     = np.array([out_HcM[15],out_HcM[16]])
epcc_LcM     = np.array([out_LcM[15],out_LcM[16]])



# -----------------------
f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(lMH_T,lM200_T,facecolor='none',edgecolors='k'  ,label='Total sample')
ax[0].scatter(lMH_L,lM200_L,facecolor='none',edgecolors='C0'  ,label='$z < 0.13$')
ax[0].scatter(lMH_H,lM200_H,facecolor='none',edgecolors='C3'  ,label='$z \geq 0.13$')

ax[0].errorbar(lMH_T,lM200_T,yerr=elM200_T,fmt = 'none',ecolor='k')
ax[0].errorbar(lMH_L,lM200_L,yerr=elM200_T,fmt = 'none',ecolor='C0')
ax[0].errorbar(lMH_H,lM200_H,yerr=elM200_T,fmt = 'none',ecolor='C3')


ax[0].plot([12.3,16],[12.3,16],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_ylabel('$\log (M_{200})$')
ax[0].axis([12.3,15.0,12.3,15.])

ax[1].scatter(lMH_T  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
ax[1].scatter(lMH_L  ,M200_L/MH_L,facecolor='none',edgecolors='C0',label='$z < 0.13$')
ax[1].scatter(lMH_H  ,M200_H/MH_H,facecolor='none',edgecolors='C3',label='$z \geq 0.13$')

ax[1].errorbar(lMH_T  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1].errorbar(lMH_L  ,M200_L/MH_L,yerr=eM200_L/MH_L,fmt = 'none',ecolor='C0')
ax[1].errorbar(lMH_H  ,M200_H/MH_H,yerr=eM200_H/MH_H,fmt = 'none',ecolor='C3')

ax[1].plot([12.3,15.5],[1.,1.],'C7--')
ax[1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1].set_ylabel(r'$M_{200}/\langle M_{AM} \rangle$')

plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/Mhalo_M200_Mbin_z.png',bbox_inches='tight')

# -----------------------
f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].plot(lMH_TcM,lM200_TcM,'ko'  ,label='Total sample')
ax[0].plot(lMH_LcM,lM200_LcM,'C0o'  ,label='$z < 0.13$')
ax[0].plot(lMH_HcM,lM200_HcM,'C3o'  ,label='$z \geq 0.13$')

ax[0].errorbar(lMH_TcM,lM200_TcM,yerr=elM200_TcM,fmt = 'none',ecolor='k')
ax[0].errorbar(lMH_LcM,lM200_LcM,yerr=elM200_TcM,fmt = 'none',ecolor='C0')
ax[0].errorbar(lMH_HcM,lM200_HcM,yerr=elM200_TcM,fmt = 'none',ecolor='C3')


ax[0].plot([12.3,16],[12.3,16],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_ylabel('$\log (M_{200})$')
ax[0].axis([12.3,15.0,12.3,15.])

ax[1].plot(lMH_TcM  ,M200_TcM/MH_TcM,'ko',label='Total sample')
ax[1].plot(lMH_LcM  ,M200_LcM/MH_LcM,'C0o',label='$z < 0.13$')
ax[1].plot(lMH_HcM  ,M200_HcM/MH_HcM,'C3o',label='$z \geq 0.13$')

ax[1].errorbar(lMH_TcM  ,M200_TcM/MH_TcM,yerr=eM200_TcM/MH_TcM,fmt = 'none',ecolor='k')
ax[1].errorbar(lMH_LcM  ,M200_LcM/MH_LcM,yerr=eM200_LcM/MH_LcM,fmt = 'none',ecolor='C0')
ax[1].errorbar(lMH_HcM  ,M200_HcM/MH_HcM,yerr=eM200_HcM/MH_HcM,fmt = 'none',ecolor='C3')

ax[1].plot([12.3,15.5],[1.,1.],'C7--')
ax[1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1].set_ylabel(r'$M_{200}/\langle M_{AM} \rangle$')

plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/Mhalo_M200_Mbin_z_cM.png',bbox_inches='tight')



# -----------------------
plt.figure()


plt.scatter(out_T[9]  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(out_L[9]  ,M200_L/MH_L,facecolor='none',edgecolors='C0',label='$z < 0.13$')
plt.scatter(out_H[9]  ,M200_H/MH_H,facecolor='none',edgecolors='C3',label='$z \geq 0.13$')

plt.errorbar(out_T[9]  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
plt.errorbar(out_L[9]  ,M200_L/MH_L,yerr=eM200_L/MH_L,fmt = 'none',ecolor='C0')
plt.errorbar(out_H[9]  ,M200_H/MH_H,yerr=eM200_H/MH_H,fmt = 'none',ecolor='C3')

plt.legend(loc= 2,frameon = False,fontsize = 13,ncol=2)


# plt.legend(frameon = False)
plt.plot([0.06,0.175],[1.,1.],'C7--')
plt.xlabel(r'$\langle z \rangle$')
plt.ylabel(r'$M_{200}/\langle M_{AM} \rangle$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/ratio_MH_z2.png',bbox_inches='tight')


# -----------------------
plt.figure()


plt.plot(out_TcM[9]  ,M200_TcM/MH_TcM,'ko',label='Total sample')
plt.plot(out_LcM[9]  ,M200_LcM/MH_LcM,'C0o',label='$z < 0.13$')
plt.plot(out_HcM[9]  ,M200_HcM/MH_HcM,'C3o',label='$z \geq 0.13$')

plt.errorbar(out_TcM[9]  ,M200_TcM/MH_TcM,yerr=eM200_TcM/MH_TcM,fmt = 'none',ecolor='k')
plt.errorbar(out_LcM[9]  ,M200_LcM/MH_LcM,yerr=eM200_LcM/MH_LcM,fmt = 'none',ecolor='C0')
plt.errorbar(out_HcM[9]  ,M200_HcM/MH_HcM,yerr=eM200_HcM/MH_HcM,fmt = 'none',ecolor='C3')

plt.legend(loc= 2,frameon = False,fontsize = 13,ncol=2)


# plt.legend(frameon = False)
plt.plot([0.06,0.175],[1.,1.],'C7--')
plt.xlabel(r'$\langle z \rangle$')
plt.ylabel(r'$M_{200}/\langle M_{AM} \rangle$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/ratio_MH_z2_cM.png',bbox_inches='tight')



# -----------------------
plt.figure()

plt.scatter(lMH_T,pcc_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lMH_L,pcc_L,facecolor='none',edgecolors='C0',label='$z < 0.13$')
plt.scatter(lMH_H,pcc_H,facecolor='none',edgecolors='C3',label='$z \geq 0.13$')

plt.errorbar(lMH_T,pcc_T,yerr=epcc_T,fmt = 'none',ecolor='k')
plt.errorbar(lMH_L,pcc_L,yerr=epcc_L,fmt = 'none',ecolor='C0')
plt.errorbar(lMH_H,pcc_H,yerr=epcc_H,fmt = 'none',ecolor='C3')

plt.plot(lMH_simu,pcc_simu,'C3')
plt.axis([12.5,14.5,0.1,1.3])
plt.legend(frameon = False,ncol = 2)

plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$p_{cc}$')
plt.savefig(r'/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/Mhalo_pcc_z.png',bbox_inches='tight')

# -----------------------
plt.figure()

plt.plot(lMH_TcM,pcc_TcM,'ko',label='Total sample')
plt.plot(lMH_LcM,pcc_LcM,'C0o',label='$z < 0.13$')
plt.plot(lMH_HcM,pcc_HcM,'C3o',label='$z \geq 0.13$')

plt.errorbar(lMH_TcM,pcc_TcM,yerr=epcc_TcM,fmt = 'none',ecolor='k')
plt.errorbar(lMH_LcM,pcc_LcM,yerr=epcc_LcM,fmt = 'none',ecolor='C0')
plt.errorbar(lMH_HcM,pcc_HcM,yerr=epcc_HcM,fmt = 'none',ecolor='C3')

plt.plot(lMH_simu,pcc_simu,'C3')
plt.axis([12.5,14.5,0.1,1.3])
plt.legend(frameon = False,ncol = 2)

plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$p_{cc}$')
plt.savefig(r'/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/Mhalo_pcc_z_cM.png',bbox_inches='tight')


# -----------------------

ratioz,eratioz     = ratio(M200_H,eM200_H,M200_L,eM200_L)
ratiozc,eratiozc     = ratio(M200_HcM,eM200_HcM,M200_LcM,eM200_LcM)

plt.figure()

plt.scatter(lMH_T  ,ratioz,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lMH_TcM  ,ratiozc,facecolor='none',edgecolors='C9',label='C-sample')

plt.errorbar(lMH_T  ,ratioz,yerr=eratioz,fmt = 'none',ecolor='k')
plt.errorbar(lMH_TcM  ,ratiozc,yerr=eratiozc,fmt = 'none',ecolor='C9')

plt.legend(loc=2,frameon = False,fontsize = 14)
# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_c_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_c_t[j]+eratio_c_t[1][j],ratio_c_t[j]-eratio_c_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,15.],[1.,1.],'C7--')
plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$M^H_{200}/M^L_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newanalysis/ratio_z.png',bbox_inches='tight')



# -----------------------
