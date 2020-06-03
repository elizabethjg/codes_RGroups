import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from profiles_fit import *

out = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin.out').T
out_zH = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_zH.out').T
out_zL = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_zL.out').T
out_mzL = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_mzL.out').T
out_mzH = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_mzH.out').T

lMH = out[4]
lMH_mzH = out_mzH[4]
lMH_mzL = out_mzL[4]
lMH_zL = out_zL[4]
lMH_zH = out_zH[4]

M200 = out[6]
M200_zL = out_zL[6]
M200_zH = out_zH[6]
M200_mzH = out_mzH[6]
M200_mzL = out_mzL[6]

eM200 = np.array([out[7],out[8]])
eM200_zH = np.array([out_zH[7],out_zH[8]])
eM200_zL = np.array([out_zL[7],out_zL[8]])
eM200_mzL = np.array([out_mzL[7],out_mzL[8]])
eM200_mzH = np.array([out_mzH[7],out_mzH[8]])

pcc = out[9]
pcc_zL = out_zL[9]
pcc_zH = out_zH[9]
pcc_mzH = out_mzH[9]
pcc_mzL = out_mzL[9]

epcc     = np.array([out[10],out[11]])
epcc_zL  = np.array([out_zL[10],out_zL[11]])
epcc_zH  = np.array([out_zH[10],out_zH[11]])
epcc_mzH = np.array([out_mzH[10],out_mzH[11]])
epcc_mzL = np.array([out_mzL[10],out_mzL[11]])

N = out[2]
mN1 = N==1
mN23 = (N>1)*(N<4)
mN410 = (N>3)*(N<11)
mN11M = N>10

# -----------------------
plt.scatter(lMH,M200,facecolor='none',edgecolors='C3')
plt.scatter(lMH_zH,M200_zH,facecolor='none',edgecolors='C1',label = 'z > m(z)')
plt.scatter(lMH_zL,M200_zL,facecolor='none',edgecolors='C0',label = 'z < m(z)')

plt.errorbar(lMH,M200,yerr=eM200,fmt = 'none',ecolor='C3')
plt.errorbar(lMH_zH,M200_zH,yerr=eM200_zH,fmt = 'none',ecolor='C1')
plt.errorbar(lMH_zL,M200_zL,yerr=eM200_zL,fmt = 'none',ecolor='C0')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------

# -----------------------
plt.scatter(lMH,M200,facecolor='none',edgecolors='C3')
plt.scatter(lMH_mzH,M200_mzH,facecolor='none',edgecolors='C1',label = 'z > m(z)')
plt.scatter(lMH_mzL,M200_mzL,facecolor='none',edgecolors='C0',label = 'z < m(z)')

plt.errorbar(lMH,M200,yerr=eM200,fmt = 'none',ecolor='C3')
plt.errorbar(lMH_mzH,M200_mzH,yerr=eM200_mzH,fmt = 'none',ecolor='C1')
plt.errorbar(lMH_mzL,M200_mzL,yerr=eM200_mzL,fmt = 'none',ecolor='C0')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------

# -----------------------
plt.figure()
plt.scatter(lMH[mN1],M200[mN1],facecolor='none',edgecolors='C6',label='N = 1')
plt.errorbar(lMH[mN1],M200[mN1],yerr=eM200[:,mN1],fmt = 'none',ecolor='C6')

plt.scatter(lMH[mN23],M200[mN23],facecolor='none',edgecolors='C7',label='N = 2-3')
plt.errorbar(lMH[mN23],M200[mN23],yerr=eM200[:,mN23],fmt = 'none',ecolor='C7')

plt.scatter(lMH[mN410],M200[mN410],facecolor='none',edgecolors='C8',label='N = 4-10')
plt.errorbar(lMH[mN410],M200[mN410],yerr=eM200[:,mN410],fmt = 'none',ecolor='C8')

plt.scatter(lMH[mN11M],M200[mN11M],facecolor='none',edgecolors='C9',label='N > 10')
plt.errorbar(lMH[mN11M],M200[mN11M],yerr=eM200[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------
##########################
# -----------------------
plt.figure()
plt.plot(lMH_mzH[mN1],M200_mzH[mN1],'C6^',label='N = 1')
plt.errorbar(lMH_mzH[mN1],M200_mzH[mN1],yerr=eM200_mzH[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(lMH_mzH[mN23],M200_mzH[mN23],'C7^',label='N = 2-3')
plt.errorbar(lMH_mzH[mN23],M200_mzH[mN23],yerr=eM200_mzH[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(lMH_mzH[mN410],M200_mzH[mN410],'C8^',label='N = 4-10')
plt.errorbar(lMH_mzH[mN410],M200_mzH[mN410],yerr=eM200_mzH[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(lMH_mzH[mN11M],M200_mzH[mN11M],'C9^',label='N > 10')
plt.errorbar(lMH_mzH[mN11M],M200_mzH[mN11M],yerr=eM200_mzH[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------

# -----------------------
plt.figure()
plt.plot(lMH_mzL[mN1],M200_mzL[mN1],'C6v',label='N = 1')
plt.errorbar(lMH_mzL[mN1],M200_mzL[mN1],yerr=eM200_mzL[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(lMH_mzL[mN23],M200_mzL[mN23],'C7v',label='N = 2-3')
plt.errorbar(lMH_mzL[mN23],M200_mzL[mN23],yerr=eM200_mzL[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(lMH_mzL[mN410],M200_mzL[mN410],'C8v',label='N = 4-10')
plt.errorbar(lMH_mzL[mN410],M200_mzL[mN410],yerr=eM200_mzL[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(lMH_mzL[mN11M],M200_mzL[mN11M],'C9v',label='N > 10')
plt.errorbar(lMH_mzL[mN11M],M200_mzL[mN11M],yerr=eM200_mzL[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------
##########################

##########################
# -----------------------
plt.figure()
plt.plot(lMH_zH[mN1],M200_zH[mN1],'C6^',label='N = 1')
plt.errorbar(lMH_zH[mN1],M200_zH[mN1],yerr=eM200_zH[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(lMH_zH[mN23],M200_zH[mN23],'C7^',label='N = 2-3')
plt.errorbar(lMH_zH[mN23],M200_zH[mN23],yerr=eM200_zH[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(lMH_zH[mN410],M200_zH[mN410],'C8^',label='N = 4-10')
plt.errorbar(lMH_zH[mN410],M200_zH[mN410],yerr=eM200_zH[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(lMH_zH[mN11M],M200_zH[mN11M],'C9^',label='N > 10')
plt.errorbar(lMH_zH[mN11M],M200_zH[mN11M],yerr=eM200_zH[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------

# -----------------------
plt.figure()
plt.plot(lMH_zL[mN1],M200_zL[mN1],'C6v',label='N = 1')
plt.errorbar(lMH_zL[mN1],M200_zL[mN1],yerr=eM200_zL[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(lMH_zL[mN23],M200_zL[mN23],'C7v',label='N = 2-3')
plt.errorbar(lMH_zL[mN23],M200_zL[mN23],yerr=eM200_zL[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(lMH_zL[mN410],M200_zL[mN410],'C8v',label='N = 4-10')
plt.errorbar(lMH_zL[mN410],M200_zL[mN410],yerr=eM200_zL[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(lMH_zL[mN11M],M200_zL[mN11M],'C9v',label='N > 10')
plt.errorbar(lMH_zL[mN11M],M200_zL[mN11M],yerr=eM200_zL[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.5,14.6],[12.5,14.6],'C7--')
plt.legend()
plt.xlabel('$\log (M_{HALO})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.5,14.5,12.5,14.5])
# -----------------------
##########################

##########################
# -----------------------
plt.figure()
plt.plot(M200_mzH[mN1],M200_mzL[mN1],'C6o',label='N = 1')
plt.errorbar(M200_mzH[mN1],M200_mzL[mN1],yerr=eM200_mzL[:,mN1],xerr=eM200_mzH[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(M200_mzH[mN23],M200_mzL[mN23],'C7o',label='N = 2-3')
plt.errorbar(M200_mzH[mN23],M200_mzL[mN23],yerr=eM200_mzL[:,mN23],xerr=eM200_mzH[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(M200_mzH[mN410],M200_mzL[mN410],'C8o',label='N = 4-10')
plt.errorbar(M200_mzH[mN410],M200_mzL[mN410],yerr=eM200_mzL[:,mN410],xerr=eM200_mzH[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(M200_mzH[mN11M],M200_mzL[mN11M],'C9o',label='N > 10')
plt.errorbar(M200_mzH[mN11M],M200_mzL[mN11M],yerr=eM200_mzL[:,mN11M],xerr=eM200_mzH[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.0,14.6],[12.0,14.6],'C7--')
plt.legend()
plt.ylabel('$\log (M^{zL}_{200})$')
plt.xlabel('$\log (M^{zH}_{200})$')
plt.axis([12.0,14.5,12.0,14.5])
# -----------------------
M200_mzH2 = 10**M200_mzH 
M200_mzL2 = 10**M200_mzL

eM200_mzH2 = [eM200_mzH[0]*M200_mzH2*np.log(10.),eM200_mzH[0]*M200_mzH2*np.log(10.)]
eM200_mzL2 = [eM200_mzL[0]*M200_mzL2*np.log(10.),eM200_mzL[0]*M200_mzL2*np.log(10.)]

ratio  = M200_mzH2/M200_mzL2
eratio = np.sqrt((eM200_mzH2/M200_mzL2)**2 + ((M200_mzH2*eM200_mzL2)/(M200_mzL2**2))**2)

plt.figure()
plt.plot(M200[mN1],ratio[mN1],'C6o',label='N = 1')
plt.errorbar(M200[mN1],ratio[mN1],yerr=eratio[:,mN1],xerr=eM200[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(M200[mN23],ratio[mN23],'C7o',label='N = 2-3')
plt.errorbar(M200[mN23],ratio[mN23],yerr=eratio[:,mN23],xerr=eM200[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(M200[mN410],ratio[mN410],'C8o',label='N = 4-10')
plt.errorbar(M200[mN410],ratio[mN410],yerr=eratio[:,mN410],xerr=eM200[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(M200[mN11M],ratio[mN11M],'C9o',label='N = 1')
plt.errorbar(M200[mN11M],ratio[mN11M],yerr=eratio[:,mN11M],xerr=eM200[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.0,14.6],[1.,1.],'C7--')
plt.legend()
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^{zH}/M^{zL}$')
# plt.axis([12.0,14.5,12.0,14.5])


# -----------------------
plt.figure()
plt.plot(M200_zH[mN1],M200_zL[mN1],'C6o',label='N = 1')
plt.errorbar(M200_zH[mN1],M200_zL[mN1],yerr=eM200_zL[:,mN1],xerr=eM200_zH[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(M200_zH[mN23],M200_zL[mN23],'C7o',label='N = 2-3')
plt.errorbar(M200_zH[mN23],M200_zL[mN23],yerr=eM200_zL[:,mN23],xerr=eM200_zH[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(M200_zH[mN410],M200_zL[mN410],'C8o',label='N = 4-10')
plt.errorbar(M200_zH[mN410],M200_zL[mN410],yerr=eM200_zL[:,mN410],xerr=eM200_zH[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(M200_zH[mN11M],M200_zL[mN11M],'C9o',label='N > 10')
plt.errorbar(M200_zH[mN11M],M200_zL[mN11M],yerr=eM200_zL[:,mN11M],xerr=eM200_zH[:,mN11M],fmt = 'none',ecolor='C9')

plt.plot([12.0,14.6],[12.0,14.6],'C7--')
plt.legend()
plt.ylabel('$\log (M^{zL}_{200})$')
plt.xlabel('$\log (M^{zH}_{200})$')
plt.axis([12.0,14.5,12.0,14.5])
# -----------------------
###################################
# -----------------------
plt.figure()
plt.plot(M200[mN1],pcc[mN1],'C6o',label='N = 1')
plt.errorbar(M200[mN1],pcc[mN1],yerr=epcc[:,mN1],xerr=eM200[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(M200[mN23],pcc[mN23],'C7o',label='N = 2-3')
plt.errorbar(M200[mN23],pcc[mN23],yerr=epcc[:,mN23],xerr=eM200[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(M200[mN410],pcc[mN410],'C8o',label='N = 4-10')
plt.errorbar(M200[mN410],pcc[mN410],yerr=epcc[:,mN410],xerr=eM200[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(M200[mN11M],pcc[mN11M],'C9o',label='N > 10')
plt.errorbar(M200[mN11M],pcc[mN11M],yerr=epcc[:,mN11M],xerr=eM200[:,mN11M],fmt = 'none',ecolor='C9')

plt.legend()
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.axis([12.0,14.5,0.4,1.0])
# -----------------------

# -----------------------
plt.figure()
plt.plot(M200_mzH[mN1],pcc_mzH[mN1],'C6^',label='N = 1')
plt.errorbar(M200_mzH[mN1],pcc_mzH[mN1],yerr=epcc_mzH[:,mN1],xerr=eM200_mzH[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(M200_mzH[mN23],pcc_mzH[mN23],'C7^',label='N = 2-3')
plt.errorbar(M200_mzH[mN23],pcc_mzH[mN23],yerr=epcc_mzH[:,mN23],xerr=eM200_mzH[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(M200_mzH[mN410],pcc_mzH[mN410],'C8^',label='N = 4-10')
plt.errorbar(M200_mzH[mN410],pcc_mzH[mN410],yerr=epcc_mzH[:,mN410],xerr=eM200_mzH[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(M200_mzH[mN11M],pcc_mzH[mN11M],'C9^',label='N > 10')
plt.errorbar(M200_mzH[mN11M],pcc_mzH[mN11M],yerr=epcc_mzH[:,mN11M],xerr=eM200_mzH[:,mN11M],fmt = 'none',ecolor='C9')

plt.legend()
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.axis([12.0,14.5,0.4,1.0])
# -----------------------

plt.figure()
plt.plot(M200_mzL[mN1],pcc_mzL[mN1],'C6v',label='N = 1')
plt.errorbar(M200_mzL[mN1],pcc_mzL[mN1],yerr=epcc_mzL[:,mN1],xerr=eM200_mzL[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(M200_mzL[mN23],pcc_mzL[mN23],'C7v',label='N = 2-3')
plt.errorbar(M200_mzL[mN23],pcc_mzL[mN23],yerr=epcc_mzL[:,mN23],xerr=eM200_mzL[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(M200_mzL[mN410],pcc_mzL[mN410],'C8v',label='N = 4-10')
plt.errorbar(M200_mzL[mN410],pcc_mzL[mN410],yerr=epcc_mzL[:,mN410],xerr=eM200_mzL[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(M200_mzL[mN11M],pcc_mzL[mN11M],'C9v',label='N > 10')
plt.errorbar(M200_mzL[mN11M],pcc_mzL[mN11M],yerr=epcc_mzL[:,mN11M],xerr=eM200_mzL[:,mN11M],fmt = 'none',ecolor='C9')

plt.legend()
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.axis([12.0,14.5,0.4,1.0])
# -----------------------

####################################
# -----------------------
plt.figure()
plt.plot(pcc_mzL[mN1],pcc_mzH[mN1],'C6^',label='N = 1')
plt.errorbar(pcc_mzL[mN1],pcc_mzH[mN1],yerr=epcc_mzH[:,mN1],xerr=epcc_mzL[:,mN1],fmt = 'none',ecolor='C6')

plt.plot(pcc_mzL[mN23],pcc_mzH[mN23],'C7^',label='N = 2-3')
plt.errorbar(pcc_mzL[mN23],pcc_mzH[mN23],yerr=epcc_mzH[:,mN23],xerr=epcc_mzL[:,mN23],fmt = 'none',ecolor='C7')

plt.plot(pcc_mzL[mN410],pcc_mzH[mN410],'C8^',label='N = 4-10')
plt.errorbar(pcc_mzL[mN410],pcc_mzH[mN410],yerr=epcc_mzH[:,mN410],xerr=epcc_mzL[:,mN410],fmt = 'none',ecolor='C8')

plt.plot(pcc_mzL[mN11M],pcc_mzH[mN11M],'C9^',label='N > 10')
plt.errorbar(pcc_mzL[mN11M],pcc_mzH[mN11M],yerr=epcc_mzH[:,mN11M],xerr=epcc_mzL[:,mN11M],fmt = 'none',ecolor='C9')

plt.legend()
plt.xlabel('$p^{zL}_{cc}$')
plt.ylabel('$p^{zH}_{cc}$')
