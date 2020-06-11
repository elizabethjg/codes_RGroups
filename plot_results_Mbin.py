import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from profiles_fit import *
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 16})
plt.rc('font', family='serif', size='16.0')

def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin.out').T
outH   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_mzH.out').T
outL   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_mzL.out').T

outc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_conM.out').T
outHc = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_mzH_conM.out').T
outLc = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles/Lens_Mbin_mzL_conM.out').T

lMH  = np.log10(out[6]*1.e14)
lMHH = np.log10(outH[6]*1.e14)
lMHL = np.log10(outL[6]*1.e14)

lMHc  = np.log10(outc[6]*1.e14)
lMHHc = np.log10(outHc[6]*1.e14)
lMHLc = np.log10(outLc[6]*1.e14)

M200  = out[11] 
M200H = outH[11]
M200L = outL[11]

M200c  = outc[11] 
M200Hc = outHc[11]
M200Lc = outLc[11]

eM200   = [out[12],out[13]]
eM200H  = [outH[12],outH[13]]
eM200L  = [outL[12],outL[13]]

eM200c  = [outc[12],outc[13]]
eM200Hc = [outHc[12],outHc[13]]
eM200Lc = [outLc[12],outLc[13]]

ratio_LH, eratio_LH = ratio(M200H,eM200H,M200L,eM200L)
ratio_LHc, eratio_LHc = ratio(M200Hc,eM200Hc,M200Lc,eM200Lc)
ratio_M, eratio_M = ratio(M200c,eM200c,M200,eM200)

lM200  = np.log10(out[11]*1.e14)
lM200H = np.log10(outH[11]*1.e14)
lM200L = np.log10(outL[11]*1.e14)

lM200c  = np.log10(outc[11]*1.e14)
lM200Hc = np.log10(outHc[11]*1.e14)
lM200Lc = np.log10(outLc[11]*1.e14)

elM200   = [out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))]
elM200H  = [outH[12]/(outH[11]*np.log(10.)),outH[13]/(outH[11]*np.log(10.))]
elM200L  = [outL[12]/(outL[11]*np.log(10.)),outL[13]/(outL[11]*np.log(10.))]

elM200c = [outc[12]/(outc[11]*np.log(10.)),outc[13]/(outc[11]*np.log(10.))]
elM200Hc = [outHc[12]/(outHc[11]*np.log(10.)),outHc[13]/(outHc[11]*np.log(10.))]
elM200Lc = [outLc[12]/(outLc[11]*np.log(10.)),outLc[13]/(outLc[11]*np.log(10.))]

pcc   = out[14]
pccH   = outH[14]
pccL   = outL[14]

pcc_c  = outc[14]
pccH_c = outHc[14]
pccL_c = outLc[14]

epcc     = np.array([out[15],out[16]])
epccH    = np.array([outH[15],outH[16]])
epccL    = np.array([outL[15],outL[16]])

epcc_c   = np.array([outc[15],outc[16]])
epccH_c   = np.array([outHc[15],outHc[16]])
epccL_c   = np.array([outLc[15],outLc[16]])

ratio_pcc, eratio_pcc = ratio(pcc_c,epcc_c,pcc,epcc)

# -----------------------
plt.figure()
plt.scatter(lMH,lM200,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lMHc,lM200c,facecolor='none',edgecolors='C9', label = u'$C-$sample')

# plt.plot(lMHH,lM200H,'k^')
# plt.plot(lMHHc,lM200Hc,'C4^')

# plt.plot(lMHL,lM200L,'kv')
# plt.plot(lMHLc,lM200Lc,'C4v')

plt.errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='k')
plt.errorbar(lMHc,lM200c,yerr=elM200c,fmt = 'none',ecolor='C9')

# plt.errorbar(lMHH,lM200H,yerr=elM200H,fmt = 'none',ecolor='k')
# plt.errorbar(lMHHc,lM200Hc,yerr=elM200Hc,fmt = 'none',ecolor='C9')

# plt.errorbar(lMHL,lM200L,yerr=elM200L,fmt = 'none',ecolor='k')
# plt.errorbar(lMHLc,lM200Lc,yerr=elM200Lc,fmt = 'none',ecolor='C9')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False,fontsize=14)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.6,12.3,14.6])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_M200_Mbin.pdf',bbox_inches='tight')
# -----------------------

plt.figure()
plt.scatter(lM200,pcc,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lM200c,pcc_c,facecolor='none',edgecolors='C9', label = u'$C-$sample')

plt.errorbar(lM200,pcc,xerr=elM200,yerr=epcc,fmt = 'none',ecolor='k')
plt.errorbar(lM200c,pcc_c,xerr=elM200c,yerr=epcc_c,fmt = 'none',ecolor='C9')


# plt.legend(frameon = False)
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/M200_pcc_Mbin.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.scatter(lM200,ratio_pcc,facecolor='none',edgecolors='k', label = u'$C-$sample')
plt.errorbar(lM200,ratio_pcc,xerr=elM200,yerr=eratio_pcc,fmt = 'none',ecolor='k')


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p^c_{cc}/p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratiopcc_Mbin.pdf',bbox_inches='tight')

# -----------------------
plt.figure()
plt.scatter(lM200,ratio_M,facecolor='none',edgecolors='k', label = u'$C-$sample')
plt.errorbar(lM200,ratio_M,xerr=elM200,yerr=eratio_M,fmt = 'none',ecolor='k')


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^c_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioc_Mbin.pdf',bbox_inches='tight')

'''
# -----------------------
plt.figure()
plt.scatter(lM200,ratio_LH,facecolor='none',edgecolors='k', label = 'Total sample')
plt.scatter(lM200c,ratio_LHc,facecolor='none',edgecolors='C9', label = u'$C-$sample')

plt.errorbar(lM200,ratio_LH,xerr=elM200,yerr=eratio_LH,fmt = 'none',ecolor='k')
plt.errorbar(lM200c,ratio_LHc,xerr=elM200c,yerr=eratio_LHc,fmt = 'none',ecolor='C9')

plt.yscale('log')
plt.legend(frameon = False)
# plt.legend(frameon = False)
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^H_{200}/M^L_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioLH_Mbin.pdf',bbox_inches='tight')


'''
