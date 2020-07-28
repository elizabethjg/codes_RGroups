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

out_T   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin.out').T
out4M   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_4MzH.out').T
outc_T  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin_cM.out').T

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Nbin.out').T
outc    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Nbin_cM.out').T

lMdyn,sdisp     = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/median_dyn').T
lMdyn_c,sdisp_c = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/median_dyn_c').T

lMH_simu = pcc_s[2]
pcc_simu = 1-(pcc_s[-1]/100.)

MH   = out[6]
MH4M = out4M[6]
MH_T = out_T[6]

lMH  = np.log10(out[6]*1.e14)
lMHc  = np.log10(outc[6]*1.e14)

lMH_T  = np.log10(out_T[6]*1.e14)
lMHc_T  = np.log10(outc_T[6]*1.e14)

M2004M  = out4M[11] 
M200  = out[11] 
M200c  = outc[11] 

M200_T  = out_T[11] 
M200c_T  = outc_T[11] 

eM2004M   = np.array([out4M[12],out4M[13]])
eM200   = np.array([out[12],out[13]])
eM200c   = np.array([outc[12],outc[13]])

eM200_T   = np.array([out_T[12],out_T[13]])
eM200c_T  = np.array([outc_T[12],outc_T[13]])

ratio_c, eratio_c = ratio(M200c,eM200c,M200,eM200)
ratio_cT, eratio_cT = ratio(M200c_T,eM200c_T,M200_T,eM200_T)

lM200_T  = np.log10(out_T[11]*1.e14)
elM200_T = np.array([out_T[12]/(out_T[11]*np.log(10.)),out_T[13]/(out_T[11]*np.log(10.))])

lM200  = np.log10(out[11]*1.e14)
lM200c  = np.log10(outc[11]*1.e14)

elM200   = np.array([out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))])
elM200c   = np.array([outc[12]/(outc[11]*np.log(10.)),outc[13]/(outc[11]*np.log(10.))])

pcc   = out[14]
pcc_T = out_T[14]
pcc_c   = outc[14]

epcc     = np.array([out[15],out[16]])
epcc_T     = np.array([out_T[15],out_T[16]])
epcc_c     = np.array([outc[15],outc[16]])

ratio_pcc,eratio_pcc     = ratio(pcc_c,epcc_c,pcc,epcc)


N = out[3]
mN1 = N==1
mN23 = (N>1)*(N<4)
mN4M = (N>3)

H        = cosmo.H(out[9][mN4M]).value/(1.0e3*pc) #H at z_pair s-1
M200_SIS  = ((2.*(sdisp*1.e3)**3)/((50**0.5)*G*H))/(Msun)

H        = cosmo.H(outc[9][mN4M]).value/(1.0e3*pc) #H at z_pair s-1
M200_SISc = ((2.*(sdisp_c*1.e3)**3)/((50**0.5)*G*H))/(Msun)

lMdyn2  = np.log10(M200_SIS)
lMdyn_c2 = np.log10(M200_SISc)



# -----------------------
plt.figure()

plt.scatter(lMdyn2,lM200[mN4M],facecolor='none',edgecolors='C8')
plt.errorbar(lMdyn2,lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8',label='$  N_{GAL} \geq 4$')

# plt.plot(lMdyn_c2,lM200c[mN4M],'C8^')
# plt.errorbar(lMdyn_c2,lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')

plt.xlabel('$\log (M_{dyn})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([13.,14.2,12.9,14.39])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/Mdyn_M200_NMbin_SIS.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
# plt.plot(100,100,'k^', label = '$C-$sample')
# plt.scatter(100,100,facecolor='none',edgecolors='k', label = 'Total sample')

plt.scatter(10.*out_T[6]    ,10.*out_T[11],facecolor='none',edgecolors='k'  ,label='Total sample')
plt.scatter(10.*out[6][mN1] ,10.*out[11][mN1],facecolor='none',edgecolors='C6'  ,label='$N_{GAL} = 1$')
plt.scatter(10.*out[6][mN23],10.*out[11][mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(10.*out[6][mN4M],10.*out[11][mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(10.*out_T[6]    ,10.*out_T[11],yerr=10.*eM200_T,fmt = 'none',ecolor='k')
plt.errorbar(10.*out[6][mN1] ,10.*out[11][mN1],yerr=10.*eM200[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(10.*out[6][mN23],10.*out[11][mN23],yerr=10.*eM200[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(10.*out[6][mN4M],10.*out[11][mN4M],yerr=10.*eM200[:,mN4M],fmt = 'none',ecolor='C8')

# plt.plot(lMHc[mN1] ,lM200c[mN1] ,'C6^')
# plt.plot(lMHc[mN23],lM200c[mN23],'C7^')
# plt.plot(lMHc[mN4M],lM200c[mN4M],'C8^')

# plt.errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6')
# plt.errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
# plt.errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([0,40.],[0,40],'C7--')
plt.legend(frameon = False,fontsize = 14)
plt.xlabel(r'$\langle M_{AM} \rangle [10^{13}M_\odot h^{-1}]$')
plt.ylabel('$M_{200} [10^{13}M_\odot h^{-1}]$')
plt.axis([0,35,0,40])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/Mhalo_M200_NMbin_lin.pdf',bbox_inches='tight')

# -----------------------


# -----------------------
plt.figure()
# plt.plot(100,100,'k^', label = '$C-$sample')
# plt.scatter(100,100,facecolor='none',edgecolors='k', label = 'Total sample')

plt.scatter(lMH_T,lM200_T,facecolor='none',edgecolors='k'  ,label='Total sample')
plt.scatter(lMH[mN1],lM200[mN1],facecolor='none',edgecolors='C6'  ,label='$N_{GAL} = 1$')
plt.scatter(lMH[mN23],lM200[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMH[mN4M],lM200[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lMH_T,lM200_T,yerr=elM200_T,fmt = 'none',ecolor='k')
plt.errorbar(lMH[mN1],lM200[mN1],yerr=elM200[:,mN1],fmt = 'none',ecolor='C6'   )
plt.errorbar(lMH[mN23],lM200[mN23],yerr=elM200[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMH[mN4M],lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8')

# plt.plot(lMHc[mN1] ,lM200c[mN1] ,'C6^')
# plt.plot(lMHc[mN23],lM200c[mN23],'C7^')
# plt.plot(lMHc[mN4M],lM200c[mN4M],'C8^')

# plt.errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6')
# plt.errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
# plt.errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False,fontsize = 14)
plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.7,12.3,14.7])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/Mhalo_M200_NMbin.pdf',bbox_inches='tight')

# -----------------------
plt.figure()

plt.scatter(lMH_T  ,ratio_cT,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lMH[mN1]  ,ratio_c[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMH[mN23] ,ratio_c[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMH[mN4M],ratio_c[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lMH_T  ,ratio_cT,yerr=eratio_cT,fmt = 'none',ecolor='k')
plt.errorbar(lMH[mN1]  ,ratio_c[mN1],yerr=eratio_c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMH[mN23] ,ratio_c[mN23],yerr=eratio_c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMH[mN4M],ratio_c[mN4M],yerr=eratio_c[:,mN4M],fmt = 'none',ecolor='C8')
plt.legend(frameon = False,fontsize = 14)
# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_c_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_c_t[j]+eratio_c_t[1][j],ratio_c_t[j]-eratio_c_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.3],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^c_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratioc_NMbin.pdf',bbox_inches='tight')

# -----------------------
plt.figure()

plt.scatter(lMH_T  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lMH[mN1]  ,(M200/MH)[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMH[mN23] ,(M200/MH)[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMH[mN4M],(M200/MH)[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lMH_T  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
plt.errorbar(lMH[mN1]  ,(M200/MH)[mN1],yerr=eM200[:,mN1]/MH[mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMH[mN23] ,(M200/MH)[mN23],yerr=eM200[:,mN23]/MH[mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMH[mN4M],(M200/MH)[mN4M],yerr=eM200[:,mN4M]/MH[mN4M],fmt = 'none',ecolor='C8')
plt.legend(frameon = False,fontsize = 14)
# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_c_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_c_t[j]+eratio_c_t[1][j],ratio_c_t[j]-eratio_c_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.5,15.],[1.,1.],'C7--')
plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel(r'$M_{200}/\langle M_{AM} \rangle$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratio_MH_NMbin.pdf',bbox_inches='tight')

# -----------------------

# -----------------------
plt.figure()

plt.scatter(out_T[9]  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(out[9][mN1]  ,(M200/MH)[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(out[9][mN23] ,(M200/MH)[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(out[9][mN4M],(M200/MH)[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')
plt.scatter(out4M[9],(M2004M/MH4M),facecolor='none',edgecolors='C8')

# plt.plot(out_T[9][lMH_T > 14.]  ,(M200_T/MH_T)[lMH_T > 14.],'C3o')
# plt.plot(out[9][lMH > 14.]  ,(M200/MH)[lMH > 14.],'C3o')

plt.errorbar(out4M[9]  ,M2004M/MH4M,yerr=eM2004M/MH4M,fmt = 'none',ecolor='C8')
plt.errorbar(out_T[9]  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
plt.errorbar(out[9][mN1]  ,(M200/MH)[mN1],yerr=eM200[:,mN1]/MH[mN1],fmt = 'none',ecolor='C6')
plt.errorbar(out[9][mN23] ,(M200/MH)[mN23],yerr=eM200[:,mN23]/MH[mN23],fmt = 'none',ecolor='C7')
plt.errorbar(out[9][mN4M],(M200/MH)[mN4M],yerr=eM200[:,mN4M]/MH[mN4M],fmt = 'none',ecolor='C8')
plt.legend(frameon = False,fontsize = 14)
# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_c_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_c_t[j]+eratio_c_t[1][j],ratio_c_t[j]-eratio_c_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([0.06,0.175],[1.,1.],'C7--')
plt.xlabel(r'$z$')
plt.ylabel(r'$M_{200}/\langle M_{AM} \rangle$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratio_MH_z.pdf',bbox_inches='tight')

# -----------------------
plt.figure()



plt.scatter(lM200_T,pcc_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lM200[mN1],pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23],pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],pcc[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200_T,pcc_T,xerr=elM200_T,yerr=epcc_T,fmt = 'none',ecolor='k')
plt.errorbar(lM200[mN1],pcc[mN1],xerr=elM200[:,mN1],yerr=epcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23],pcc[mN23],xerr=elM200[:,mN23],yerr=epcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],pcc[mN4M],xerr=elM200[:,mN4M],yerr=epcc[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot(lMH_simu,pcc_simu,'C3')

# plt.plot(lM200c[mN1] ,pcc_c[mN1]  ,'C6^',label='$N_{GAL} = 1$')
# plt.plot(lM200c[mN23],pcc_c[mN23],'C7^',label='$ 2 \leq N_{GAL} \leq 3$')
# plt.plot(lM200c[mN4M],pcc_c[mN4M],'C8^',label='$ 4 \leq N_{GAL} \leq 10$')

# plt.errorbar(lM200c[mN1] ,pcc_c[mN1] ,xerr=elM200c[:,mN1] ,yerr=epcc_c[:,mN1],fmt = 'none',ecolor='C6')
# plt.errorbar(lM200c[mN23],pcc_c[mN23],xerr=elM200c[:,mN23],yerr=epcc_c[:,mN23],fmt = 'none',ecolor='C7')
# plt.errorbar(lM200c[mN4M],pcc_c[mN4M],xerr=elM200c[:,mN4M],yerr=epcc_c[:,mN4M],fmt = 'none',ecolor='C8')


# plt.legend(frameon = False)
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/M200_pcc_NMbin.pdf',bbox_inches='tight')

# -----------------------

'''



plt.figure()

plt.scatter(lM200[mN1]  ,ratio_pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23] ,ratio_pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],ratio_pcc[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200[mN1]  ,ratio_pcc[mN1],xerr=elM200[:,mN1],yerr=eratio_pcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23] ,ratio_pcc[mN23],xerr=elM200[:,mN23],yerr=eratio_pcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],ratio_pcc[mN4M],xerr=elM200[:,mN4M],yerr=eratio_pcc[:,mN4M],fmt = 'none',ecolor='C8')

# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_pcc_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_pcc_t[j]+eratio_pcc_t[1][j],ratio_pcc_t[j]-eratio_pcc_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.3],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p^c_{cc}/p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratiopcc_NMbin.pdf',bbox_inches='tight')

# -----------------------





# -----------------------
plt.figure()
plt.plot(lMHc_t,lM200c_t,'C4')
plt.scatter(lMHc[mN1],lM200c[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMHc[mN23],lM200c[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMHc[mN4M],lM200c[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.fill_between(lMHc_t,lM200c_t-elM200c_t[0],lM200_t+elM200c_t[1],color='C4',alpha=0.5)
# plt.errorbar(lMHc_t,lM200c_t,yerr=elM200c_t,fmt = 'none',ecolor='C0')
plt.errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.2,12.3,14.2])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/Mhalo_M200_NMbin_conM.pdf',bbox_inches='tight')
# -----------------------

# -----------------------
plt.figure()

plt.scatter(lM200[mN1]  ,ratio_LH[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23] ,ratio_LH[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],ratio_LH[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200[mN1]  ,ratio_LH[mN1],xerr=elM200[:,mN1],yerr=eratio_LH[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23] ,ratio_LH[mN23],xerr=elM200[:,mN23],yerr=eratio_LH[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],ratio_LH[mN4M],xerr=elM200[:,mN4M],yerr=eratio_LH[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot(lM200c[mN1]  ,ratio_LHc[mN1],'C6^',label='$N_{GAL} = 1$')
plt.plot(lM200c[mN23] ,ratio_LHc[mN23],'C7^',label='$ 2 \leq N_{GAL} \leq 3$')
plt.plot(lM200c[mN4M],ratio_LHc[mN4M],'C8^',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200c[mN1]  ,ratio_LHc[mN1],xerr=elM200c[:,mN1],yerr=eratio_LHc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200c[mN23] ,ratio_LHc[mN23],xerr=elM200c[:,mN23],yerr=eratio_LHc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200c[mN4M],ratio_LHc[mN4M],xerr=elM200c[:,mN4M],yerr=eratio_LHc[:,mN4M],fmt = 'none',ecolor='C8')

j = np.argsort(lM200c_t)
plt.plot(lM200c_t[j],ratio_LHc_t[j],'C4')
plt.fill_between(lM200c_t[j],ratio_LHc_t[j]+eratio_LHc_t[1][j],ratio_LHc_t[j]-eratio_LHc_t[0][j],color='C4',alpha=0.5)

# plt.plot(lM200_t,ratio_LH_t,'k')
# plt.fill_between(lM200_t,ratio_LH_t-eratio_LH_t[0],ratio_LH_t+eratio_LH_t[1],color='k',alpha=0.5)


plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.3],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^H_{200}/M^L_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratioLH_NMbin.pdf',bbox_inches='tight')

# -----------------------
# -----------------------
plt.figure()

plt.scatter(lMdyn,lM200[mN4M],facecolor='none',edgecolors='C8')
plt.errorbar(lMdyn,lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8',label='$  N_{GAL} \geq 4$')

plt.plot(lMdyn_c,lM200c[mN4M],'C8^')
plt.errorbar(lMdyn_c,lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')

plt.xlabel('$\log (M_{dyn})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([13.,14.5,12.7,14.5])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/Mdyn_M200_NMbin.pdf',bbox_inches='tight')


'''
