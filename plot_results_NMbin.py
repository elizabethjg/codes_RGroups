import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from profiles_fit import *

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin.out').T
outH   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_mzH.out').T
outL   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_mzL.out').T

lMH  = np.log10(out[6]*1.e14)
lMHH = np.log10(outH[6]*1.e14)
lMHL = np.log10(outL[6]*1.e14)

M200  = out[11] 
M200H = outH[11]
M200L = outL[11]

eM200   = [out[12],out[13]]
eM200H  = [outH[12],outH[13]]
eM200L  = [outL[12],outL[13]]

ratio_LH  = M200H/M200L
eratio_LH = np.sqrt((eM200H/M200L)**2 + ((eM200L*M200H)/(M200L**2))**2)


lM200  = np.log10(out[11]*1.e14)
lM200H = np.log10(outH[11]*1.e14)
lM200L = np.log10(outL[11]*1.e14)

elM200   = np.array([out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))])
elM200H  = np.array([outH[12]/(outH[11]*np.log(10.)),outH[13]/(outH[11]*np.log(10.))])
elM200L  = np.array([outL[12]/(outL[11]*np.log(10.)),outL[13]/(outL[11]*np.log(10.))])

pcc   = out[14]
pccH   = outH[14]
pccL   = outL[14]

epcc     = np.array([out[15],out[16]])
epccH    = np.array([outH[15],outH[16]])
epccL    = np.array([outL[15],outL[16]])

N = out[3]
mN1 = N==1
mN23 = (N>1)*(N<4)
mN410 = (N>3)*(N<11)
mN11M = N>10


# -----------------------
plt.figure()
plt.scatter(lMH[mN1],lM200[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMH[mN23],lM200[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMH[mN410],lM200[mN410],facecolor='none',edgecolors='C8',label='$ 4 \leq N_{GAL} \leq 10$')
plt.scatter(lMH[mN11M],lM200[mN11M],facecolor='none',edgecolors='C9',label='$N_{GAL} \geq 11$')

plt.errorbar(lMH[mN1],lM200[mN1],yerr=elM200[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMH[mN23],lM200[mN23],yerr=elM200[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMH[mN410],lM200[mN410],yerr=elM200[:,mN410],fmt = 'none',ecolor='C8')
plt.errorbar(lMH[mN11M],lM200[mN11M],yerr=elM200[:,mN11M],fmt = 'none',ecolor='C9')



plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.6,12.3,14.6])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_M200_NMbin.pdf',bbox_inches='tight')
# -----------------------

plt.figure()

plt.scatter(lM200[mN1],pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23],pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN410],pcc[mN410],facecolor='none',edgecolors='C8',label='$ 4 \leq N_{GAL} \leq 10$')
plt.scatter(lM200[mN11M],pcc[mN11M],facecolor='none',edgecolors='C9',label='$N_{GAL} \geq 10$')

plt.errorbar(lM200[mN1],pcc[mN1],xerr=elM200[:,mN1],yerr=epcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23],pcc[mN23],xerr=elM200[:,mN23],yerr=epcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN410],pcc[mN410],xerr=elM200[:,mN410],yerr=epcc[:,mN410],fmt = 'none',ecolor='C8')
plt.errorbar(lM200[mN11M],pcc[mN11M],xerr=elM200[:,mN11M],yerr=epcc[:,mN11M],fmt = 'none',ecolor='C9')


# plt.legend(frameon = False)
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/M200_pcc_NMbin.pdf',bbox_inches='tight')

# -----------------------
plt.figure()

plt.scatter(lM200[mN1]  ,ratio_LH[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23] ,ratio_LH[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN410],ratio_LH[mN410],facecolor='none',edgecolors='C8',label='$ 4 \leq N_{GAL} \leq 10$')
plt.scatter(lM200[mN11M],ratio_LH[mN11M],facecolor='none',edgecolors='C9',label='$N_{GAL} \geq 10$')

plt.errorbar(lM200[mN1]  ,ratio_LH[mN1],xerr=elM200[:,mN1],yerr=eratio_LH[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23] ,ratio_LH[mN23],xerr=elM200[:,mN23],yerr=eratio_LH[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN410],ratio_LH[mN410],xerr=elM200[:,mN410],yerr=eratio_LH[:,mN410],fmt = 'none',ecolor='C8')
plt.errorbar(lM200[mN11M],ratio_LH[mN11M],xerr=elM200[:,mN11M],yerr=eratio_LH[:,mN11M],fmt = 'none',ecolor='C9')

plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^H_{200}/M^L_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioLH_NMbin.pdf',bbox_inches='tight')

# -----------------------

