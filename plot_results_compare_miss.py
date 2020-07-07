import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from profiles_fit import *

def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin.out').T
outc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin_cM.out').T

outy    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin_ymiss.out').T
outyc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin_cM_ymiss.out').T

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

ratio_M, eratio_M    = ratio(M200y,eM200y,M200,eM200)
ratio_Mc, eratio_Mc  = ratio(M200yc,eM200yc,M200c,eM200c)

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

ratio_p, eratio_p  = ratio(pccy,epcc_y,pcc,epcc)
ratio_pc, eratio_pc  = ratio(pcc_yc,epcc_yc,pcc,epcc)

# -----------------------
plt.figure()
plt.scatter(lM200,ratio_M,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lM200c,ratio_Mc,facecolor='none',edgecolors='C9', label = u'$C-$sample')

plt.errorbar(lM200,ratio_M,xerr=elM200,yerr=eratio_M,fmt = 'none',ecolor='k')
plt.errorbar(lM200c,ratio_Mc,xerr=elM200c,yerr=eratio_Mc,fmt = 'none',ecolor='C9')

# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^y_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratio_Mbin_centreY.pdf',bbox_inches='tight')


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
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratiopcc_Mbin_centreY.pdf',bbox_inches='tight')

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

plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/Mhalo_Mbin_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()

plt.plot(lMH,chi2,'ko')
plt.plot(lMHy,chi2y,'kx')
plt.plot(lMHc,chi2c,'C9o')
plt.plot(lMHyc,chi2yc,'C9x')
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$Chi^2$')
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/chi2_Mbin_centreY.pdf',bbox_inches='tight')

# plt.errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='k')
# plt.errorbar(lMHc,lM200c,yerr=elM200c,fmt = 'none',ecolor='C9')

# plt.errorbar(lMHy,lM200y,yerr=elM200y,fmt = 'none',ecolor='k')
# plt.errorbar(lMHyc,lM200yc,yerr=elM200yc,fmt = 'none',ecolor='C9')

# plt.plot([12.3,14.6],[12.3,14.6],'C7--')
# plt.legend(frameon = False)
# plt.xlabel('$\log (M_{AM})$')
# plt.ylabel('$\log (M_{200})$')
# plt.axis([12.3,14.6,12.3,14.6])
