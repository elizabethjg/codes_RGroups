import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from profiles_fit import *

def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin.out').T
outc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_cM.out').T

outy    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_ymiss.out').T
outyc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_cM_ymiss.out').T

chi2 = out[-1]
chi2c = outc[-1]
chi2y = outy[-1]
chi2yc = outyc[-1]


lMH  = np.log10(out[6]*1.e14)
lMHc  = np.log10(outc[6]*1.e14)

lMHy   = np.log10(outy[6]*1.e14)
lMHyc  = np.log10(outyc[6]*1.e14)

M200  = out[11] 
M200c  = outc[11] 

M200y  = outy[11] 
M200yc  = outyc[11] 

eM200   = np.array([out[12],out[13]])
eM200c  = np.array([outc[12],outc[13]])

eM200y   = np.array([outy[12],outy[13]])
eM200yc  = np.array([outyc[12],outyc[13]])

ratio_M, eratio_M    = ratio(M200y,eM200y,M200,eM200)
ratio_Mc, eratio_Mc  = ratio(M200yc,eM200yc,M200c,eM200c)

lM200    = np.log10(out[11]*1.e14)
lM200c   = np.log10(outc[11]*1.e14)
lM200y   = np.log10(outy[11]*1.e14)
lM200yc  = np.log10(outyc[11]*1.e14)

elM200   = np.array([out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))])
elM200c = np.array([outc[12]/(outc[11]*np.log(10.)),outc[13]/(outc[11]*np.log(10.))])

elM200y   = np.array([outy[12]/(outy[11]*np.log(10.)),outy[13]/(outy[11]*np.log(10.))])
elM200yc  = np.array([outyc[12]/(outyc[11]*np.log(10.)),outyc[13]/(outyc[11]*np.log(10.))])

pcc   = out[14]
pcc_c  = outc[14]
pccy   = outy[14]
pcc_yc  = outyc[14]

epcc_y     = np.array([outy[15],outy[16]])
epcc_yc    = np.array([outyc[15],outyc[16]])
epcc       = np.array([out[15],out[16]])
epcc_c     = np.array([outc[15],outc[16]])

ratio_p, eratio_p  = ratio(pccy,epcc_y,pcc,epcc)
ratio_pc, eratio_pc  = ratio(pcc_yc,epcc_yc,pcc,epcc)

N = out[3]
mN1 = N==1
mN23 = (N>1)*(N<4)
mN4M = (N>3)

# -----------------------
plt.figure()
plt.scatter(lM200[mN1],ratio_M[mN1],facecolor='none',edgecolors='C6', label = '$N_{GAL} = 1$')
plt.scatter(lM200[mN23],ratio_M[mN23],facecolor='none',edgecolors='C7', label = '$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],ratio_M[mN4M],facecolor='none',edgecolors='C8', label = '$  N_{GAL} \geq 4$')

plt.errorbar(lM200[mN1], ratio_M[mN1], xerr=elM200[:,mN1] ,yerr=eratio_M[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23],ratio_M[mN23],xerr=elM200[:,mN23],yerr=eratio_M[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],ratio_M[mN4M],xerr=elM200[:,mN4M],yerr=eratio_M[:,mN4M],fmt = 'none',ecolor='C8')

# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^y_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratio_NMbin_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.scatter(lM200c[mN1],ratio_Mc[mN1],facecolor='none',edgecolors='C6', label = '$N_{GAL} = 1$')
plt.scatter(lM200c[mN23],ratio_Mc[mN23],facecolor='none',edgecolors='C7', label = '$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200c[mN4M],ratio_Mc[mN4M],facecolor='none',edgecolors='C8', label = '$  N_{GAL} \geq 4$')

plt.errorbar(lM200c[mN1] ,ratio_Mc[mN1], xerr=elM200c[:,mN1] ,yerr=eratio_Mc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200c[mN23],ratio_Mc[mN23],xerr=elM200c[:,mN23],yerr=eratio_Mc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200c[mN4M],ratio_Mc[mN4M],xerr=elM200c[:,mN4M],yerr=eratio_Mc[:,mN4M],fmt = 'none',ecolor='C8')

# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^y_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratio_NMbin_cM_centreY.pdf',bbox_inches='tight')


# -----------------------
plt.figure()
plt.scatter(lM200,ratio_p,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lM200c,ratio_pc,facecolor='none',edgecolors='C9', label = u'$C-$sample')

plt.errorbar(lM200,ratio_p,xerr=elM200,yerr=eratio_p,fmt = 'none',ecolor='k')
plt.errorbar(lM200c,ratio_pc,xerr=elM200c,yerr=eratio_pc,fmt = 'none',ecolor='C9')

# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$pcc^y/pcc$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratiopcc_NMbin_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.scatter(lMH,lM200,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lMHc,lM200c,facecolor='none',edgecolors='C9', label = u'$C-$sample')

plt.plot(lMHy,lM200y,'ko')
plt.plot(lMHyc,lM200yc,'C9o')

plt.errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='k')
plt.errorbar(lMHc,lM200c,yerr=elM200c,fmt = 'none',ecolor='C9')

plt.errorbar(lMHy,lM200y,yerr=elM200y,fmt = 'none',ecolor='k')
plt.errorbar(lMHyc,lM200yc,yerr=elM200yc,fmt = 'none',ecolor='C9')

plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.6,12.3,14.6])

plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_NMbin_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()

plt.plot(lMH,chi2,'ko')
plt.plot(lMHy,chi2y,'kx')
plt.plot(lMHc,chi2c,'C9o')
plt.plot(lMHyc,chi2yc,'C9x')
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$Chi^2$')
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/chi2_NMbin_centreY.pdf',bbox_inches='tight')

# plt.errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='k')
# plt.errorbar(lMHc,lM200c,yerr=elM200c,fmt = 'none',ecolor='C9')

# plt.errorbar(lMHy,lM200y,yerr=elM200y,fmt = 'none',ecolor='k')
# plt.errorbar(lMHyc,lM200yc,yerr=elM200yc,fmt = 'none',ecolor='C9')

# plt.plot([12.3,14.6],[12.3,14.6],'C7--')
# plt.legend(frameon = False)
# plt.xlabel('$\log (M_{AM})$')
# plt.ylabel('$\log (M_{200})$')
# plt.axis([12.3,14.6,12.3,14.6])
