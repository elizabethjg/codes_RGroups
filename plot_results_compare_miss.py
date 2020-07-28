import sys
import numpy as np
from matplotlib import *
from astropy.io import fits

def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin.out').T
outy    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/Lens_Mbin_ymiss.out').T
outy2    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_new/old_ymiss/Lens_Mbin_ymiss.out').T


chi2 = out[-1]
chi2y = outy[-2]
chi2y2 = outy2[-1]

lMH  = np.log10(out[6]*1.e14)
lMHy   = np.log10(outy[6]*1.e14)
lMHy2  = np.log10(outy2[6]*1.e14)

M200  = out[11] 
M200y  = outy[11] 
M200y2  = outy2[11] 

eM200   = [out[12],out[13]]
eM200y   = [outy[12],outy[13]]
eM200y2  = [outy2[12],outy2[13]]

ratio_M, eratio_M    = ratio(M200y,eM200y,M200,eM200)
ratio_M2, eratio_M2    = ratio(M200y2,eM200y2,M200y,eM200y)

lM200    = np.log10(out[11]*1.e14)
lM200y   = np.log10(outy[11]*1.e14)
lM200y2  = np.log10(outy2[11]*1.e14)

elM200   = [out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))]
elM200y   = [outy[12]/(outy[11]*np.log(10.)),outy[13]/(outy[11]*np.log(10.))]
elM200y2 = [outy2[12]/(outy2[11]*np.log(10.)),outy2[13]/(outy2[11]*np.log(10.))]

pcc   = out[14]
pccy   = outy[14]
pccy2  = outy2[14]

epcc_y     = np.array([outy[15],outy[16]])
epcc_y2    = np.array([outy2[15],outy2[16]])
epcc       = np.array([out[15],out[16]])


ratio_p, eratio_p  = ratio(pccy,epcc_y,pcc,epcc)


# -----------------------
plt.figure()
plt.plot(chi2,chi2y,'o')
plt.plot([0.8,2.75],[0.8,2.75],'k--')
plt.xlabel('$\chi^2$')
plt.ylabel('$\chi^2_{Y}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/chi2_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.plot(epcc[0],epcc_y[0],'o')
plt.plot(epcc[1],epcc_y[1],'o')
plt.plot([0.11,0.25],[0.11,0.25],'k--')
plt.ylabel('$\epsilon_{Y}$')
plt.xlabel('$\epsilon$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/epcc_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.plot(elM200[0]/lM200,elM200y[0]/lM200y,'o')
plt.plot(elM200[1]/lM200,elM200y[1]/lM200y,'o')
plt.plot([0.0,0.01],[0.0,0.01],'k--')
plt.ylabel('$\epsilon_{Y}$')
plt.xlabel('$\epsilon$')
plt.axis([0.055,0.12,0.07,0.135])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/elM200_centreY.pdf',bbox_inches='tight')


# -----------------------
plt.figure()
plt.scatter(lM200,ratio_M,facecolor='none',edgecolors='k', label = 'Total sample')
plt.errorbar(lM200,ratio_M,xerr=elM200,yerr=eratio_M,fmt = 'none',ecolor='k')

# plt.scatter(lM200,ratio_M2,facecolor='none',edgecolors='C1', label = 'Total sample')
# plt.errorbar(lM200,ratio_M2,xerr=elM200,yerr=eratio_M2,fmt = 'none',ecolor='k')

# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^y_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratio_Mbin_centreY.pdf',bbox_inches='tight')


# -----------------------
plt.figure()
plt.scatter(lM200,ratio_p,facecolor='none',edgecolors='k', label = 'Total sample')
plt.errorbar(lM200,ratio_p,xerr=elM200,yerr=eratio_p,fmt = 'none',ecolor='k')


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$pcc^y/pcc$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots_newzbin/ratiopcc_Mbin_centreY.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.scatter(lMH,lM200,facecolor='none',edgecolors='k', label = 'Total sample')
plt.plot(lMHy,lM200y,'ko')
plt.errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='k')
plt.errorbar(lMHy,lM200y,yerr=elM200y,fmt = 'none',ecolor='k')


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
