import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from profiles_fit import *

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin.out').T
outc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_conM.out').T

outy    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_Yan.out').T
outyc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_conM_Yan.out').T

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

eM200   = [out[12],out[13]]
eM200c  = [outc[12],outc[13]]

eM200y   = [outy[12],outy[13]]
eM200yc  = [outyc[12],outyc[13]]

ratio  = M200/M200y
eratio = np.sqrt((eM200/M200y)**2 + ((eM200y*M200)/(M200y**2))**2)

ratio_c  = M200c/M200yc
eratio_c = np.sqrt((eM200c/M200yc)**2 + ((eM200yc*M200c)/(M200yc**2))**2)


lM200    = np.log10(out[11]*1.e14)
lM200c   = np.log10(outc[11]*1.e14)
lM200y   = np.log10(outy[11]*1.e14)
lM200yc  = np.log10(outyc[11]*1.e14)

elM200   = [out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))]
elM200c = [outc[12]/(outc[11]*np.log(10.)),outc[13]/(outc[11]*np.log(10.))]

elM200y   = [outy[12]/(outy[11]*np.log(10.)),outy[13]/(outy[11]*np.log(10.))]
elM200yc = [outyc[12]/(outyc[11]*np.log(10.)),outyc[13]/(outyc[11]*np.log(10.))]


pcc   = out[14]
pcc_c  = outc[14]
pccy   = outy[14]
pcc_yc  = outyc[14]

epcc_y     = np.array([outy[15],outy[16]])
epcc_yc    = np.array([outyc[15],outyc[16]])
epcc       = np.array([out[15],out[16]])
epcc_c     = np.array([outc[15],outc[16]])

ratio_p  = pcc/pccy
eratio_p = np.sqrt((epcc/pccy)**2 + ((epcc_y*pcc)/(pccy**2))**2)

ratio_pc  = pcc_c/pcc_yc
eratio_pc = np.sqrt((epcc_c/pcc_yc)**2 + ((epcc_yc*pcc_c)/(pcc_yc**2))**2)


# -----------------------
plt.figure()
plt.scatter(lM200,ratio,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lM200c,ratio_c,facecolor='none',edgecolors='C4', label = u'$C-$sample')

plt.errorbar(lM200,ratio,xerr=elM200,yerr=eratio,fmt = 'none',ecolor='k')
plt.errorbar(lM200c,ratio_c,xerr=elM200c,yerr=eratio_c,fmt = 'none',ecolor='C4')

plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M_{200}/M^y_{200}$')
# plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioLH_Mbin.pdf',bbox_inches='tight')


# -----------------------
plt.figure()
plt.scatter(lM200,ratio_p,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lM200c,ratio_pc,facecolor='none',edgecolors='C4', label = u'$C-$sample')

plt.errorbar(lM200,ratio_p,xerr=elM200,yerr=eratio_p,fmt = 'none',ecolor='k')
plt.errorbar(lM200c,ratio_pc,xerr=elM200c,yerr=eratio_pc,fmt = 'none',ecolor='C4')

plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$pcc/pcc^y$')
# plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioLH_Mbin.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.scatter(lMH,lM200,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lMHc,lM200c,facecolor='none',edgecolors='C4', label = u'$C-$sample')

plt.plot(lMHy,lM200y,'ko')
plt.plot(lMHyc,lM200yc,'C4o')

plt.errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='k')
plt.errorbar(lMHc,lM200c,yerr=elM200c,fmt = 'none',ecolor='C4')

plt.errorbar(lMHy,lM200y,yerr=elM200y,fmt = 'none',ecolor='k')
plt.errorbar(lMHyc,lM200yc,yerr=elM200yc,fmt = 'none',ecolor='C4')

plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.6,12.3,14.6])
