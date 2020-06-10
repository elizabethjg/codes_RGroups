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
outH   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_mzH.out').T
outL   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_mzL.out').T

outc    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_cM.out').T
outHc   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_cM_mzH.out').T
outLc   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_NMbin_cM_mzL.out').T

out_t    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_nMbin_cM.out').T
outH_t   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_nMbin_mzH.out').T
outL_t   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_nMbin_mzL.out').T

outc_t    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_nMbin_cM.out').T
outHc_t   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_nMbin_cM_mzH.out').T
outLc_t   = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_NMbin/Lens_nMbin_cM_mzL.out').T


lMH  = np.log10(out[6]*1.e14)
lMHH = np.log10(outH[6]*1.e14)
lMHL = np.log10(outL[6]*1.e14)

lMHc  = np.log10(outc[6]*1.e14)
lMHHc = np.log10(outHc[6]*1.e14)
lMHLc = np.log10(outLc[6]*1.e14)

lMH_t  = np.log10(out_t[6]*1.e14)
lMHH_t = np.log10(outH_t[6]*1.e14)
lMHL_t = np.log10(outL_t[6]*1.e14)

lMHc_t  = np.log10(outc_t[6]*1.e14)
lMHHc_t = np.log10(outHc_t[6]*1.e14)
lMHLc_t = np.log10(outLc_t[6]*1.e14)


M200  = out[11] 
M200H = outH[11]
M200L = outL[11]

M200c  = outc[11] 
M200Hc = outHc[11]
M200Lc = outLc[11]

M200_t  = out_t[11] 
M200H_t = outH_t[11]
M200L_t = outL_t[11]

M200c_t  = outc_t[11] 
M200Hc_t = outHc_t[11]
M200Lc_t = outLc_t[11]


eM200   = [out[12],out[13]]
eM200H  = [outH[12],outH[13]]
eM200L  = [outL[12],outL[13]]

eM200c   = [outc[12],outc[13]]
eM200Hc  = [outHc[12],outHc[13]]
eM200Lc  = [outLc[12],outLc[13]]

eM200_t   = [out_t[12],out_t[13]]
eM200H_t  = [outH_t[12],outH_t[13]]
eM200L_t  = [outL_t[12],outL_t[13]]

eM200c_t   = [outc_t[12],outc_t[13]]
eM200Hc_t  = [outHc_t[12],outHc_t[13]]
eM200Lc_t  = [outLc_t[12],outLc_t[13]]

ratio_LH, eratio_LH   = ratio(M200H,eM200H,M200L,eM200L)
ratio_LHc, eratio_LHc = ratio(M200Hc,eM200Hc,M200Lc,eM200Lc)
ratio_c, eratio_c = ratio(M200c,eM200c,M200,eM200)

ratio_LH_t, eratio_LH_t   = ratio(M200H_t,eM200H_t,M200L_t,eM200L_t)
ratio_LHc_t, eratio_LHc_t = ratio(M200Hc_t,eM200Hc,M200Lc_t,eM200Lc_t)
ratio_c_t, eratio_c_t     = ratio(M200c_t,eM200c_t,M200_t,eM200_t)


lM200  = np.log10(out[11]*1.e14)
lM200H = np.log10(outH[11]*1.e14)
lM200L = np.log10(outL[11]*1.e14)

lM200c  = np.log10(outc[11]*1.e14)
lM200Hc = np.log10(outHc[11]*1.e14)
lM200Lc = np.log10(outLc[11]*1.e14)

lM200_t  = np.log10(out_t[11]*1.e14)
lM200H_t = np.log10(outH_t[11]*1.e14)
lM200L_t = np.log10(outL_t[11]*1.e14)

lM200c_t  = np.log10(outc_t[11]*1.e14)
lM200Hc_t = np.log10(outHc_t[11]*1.e14)
lM200Lc_t = np.log10(outLc_t[11]*1.e14)

elM200   = np.array([out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))])
elM200H  = np.array([outH[12]/(outH[11]*np.log(10.)),outH[13]/(outH[11]*np.log(10.))])
elM200L  = np.array([outL[12]/(outL[11]*np.log(10.)),outL[13]/(outL[11]*np.log(10.))])

elM200c   = np.array([outc[12]/(outc[11]*np.log(10.)),outc[13]/(outc[11]*np.log(10.))])
elM200Hc  = np.array([outHc[12]/(outHc[11]*np.log(10.)),outHc[13]/(outHc[11]*np.log(10.))])
elM200Lc  = np.array([outLc[12]/(outLc[11]*np.log(10.)),outLc[13]/(outLc[11]*np.log(10.))])

elM200_t   = np.array([out_t[12]/(out_t[11]*np.log(10.)),out_t[13]/(out_t[11]*np.log(10.))])
elM200H_t  = np.array([outH_t[12]/(outH_t[11]*np.log(10.)),outH_t[13]/(outH_t[11]*np.log(10.))])
elM200L_t  = np.array([outL_t[12]/(outL_t[11]*np.log(10.)),outL_t[13]/(outL_t[11]*np.log(10.))])

elM200c_t  = np.array([outc_t[12]/(outc_t[11]*np.log(10.)),outc_t[13]/(outc_t[11]*np.log(10.))])
elM200Hc_t = np.array([outHc_t[12]/(outHc_t[11]*np.log(10.)),outHc_t[13]/(outHc_t[11]*np.log(10.))])
elM200Lc_t = np.array([outLc_t[12]/(outLc_t[11]*np.log(10.)),outLc_t[13]/(outLc_t[11]*np.log(10.))])

pcc   = out[14]
pccH   = outH[14]
pccL   = outL[14]

pcc_c   = outc[14]
pccH_c   = outHc[14]
pccL_c   = outLc[14]

pcc_t   = out_t[14]
pccH_t   = outH_t[14]
pccL_t   = outL_t[14]

pcc_ct   = outc_t[14]
pccH_ct   = outHc_t[14]
pccL_ct   = outLc_t[14]

epcc     = np.array([out[15],out[16]])
epccH    = np.array([outH[15],outH[16]])
epccL    = np.array([outL[15],outL[16]])

epcc_c     = np.array([outc[15],outc[16]])
epccH_c    = np.array([outHc[15],outHc[16]])
epccL_c    = np.array([outLc[15],outLc[16]])

epcc_t     = np.array([out_t[15],out_t[16]])
epccH_t    = np.array([outH_t[15],outH_t[16]])
epccL_t    = np.array([outL_t[15],outL_t[16]])

epcc_ct     = np.array([outc_t[15],outc_t[16]])
epccH_ct    = np.array([outHc_t[15],outHc_t[16]])
epccL_ct    = np.array([outLc_t[15],outLc_t[16]])

ratio_pcc,eratio_pcc     = ratio(pcc_c,epcc_c,pcc,epcc)
ratio_pcc_t,eratio_pcc_t = ratio(pcc_ct,epcc_ct,pcc_t,epcc_t)

N = out[3]
mN1 = N==1
mN23 = (N>1)*(N<4)
mN4M = (N>3)

# -----------------------
plt.figure()
plt.plot(100,100,'k^', label = '$C-$sample')
plt.scatter(100,100,facecolor='none',edgecolors='k', label = 'Total sample')
# plt.plot(lMH_t,lM200_t,'k')
plt.scatter(lMH[mN1],lM200[mN1],facecolor='none',edgecolors='C6'  )
plt.scatter(lMH[mN23],lM200[mN23],facecolor='none',edgecolors='C7')
plt.scatter(lMH[mN4M],lM200[mN4M],facecolor='none',edgecolors='C8')

# plt.fill_between(lMH_t,lM200_t-elM200_t[0],lM200_t+elM200_t[1],color='k',alpha=0.3)
# plt.errorbar(lMH_t,lM200_t,yerr=elM200_t,fmt = 'none',ecolor='C0')
plt.errorbar(lMH[mN1],lM200[mN1],yerr=elM200[:,mN1],fmt = 'none',ecolor='C6'   ,label='$N_{GAL} = 1$')
plt.errorbar(lMH[mN23],lM200[mN23],yerr=elM200[:,mN23],fmt = 'none',ecolor='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.errorbar(lMH[mN4M],lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8',label='$  N_{GAL} \geq 4$')

plt.plot(lMHc[mN1] ,lM200c[mN1] ,'C6^')
plt.plot(lMHc[mN23],lM200c[mN23],'C7^')
plt.plot(lMHc[mN4M],lM200c[mN4M],'C8^')

plt.errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.3,14.2,12.3,14.2])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_M200_NMbin.pdf',bbox_inches='tight')
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
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_M200_NMbin_conM.pdf',bbox_inches='tight')
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
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioLH_NMbin.pdf',bbox_inches='tight')

# -----------------------

# -----------------------
plt.figure()

plt.scatter(lM200[mN1]  ,ratio_c[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23] ,ratio_c[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],ratio_c[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200[mN1]  ,ratio_c[mN1],xerr=elM200[:,mN1],yerr=eratio_c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23] ,ratio_c[mN23],xerr=elM200[:,mN23],yerr=eratio_c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],ratio_c[mN4M],xerr=elM200[:,mN4M],yerr=eratio_c[:,mN4M],fmt = 'none',ecolor='C8')

# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_c_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_c_t[j]+eratio_c_t[1][j],ratio_c_t[j]-eratio_c_t[0][j],color='C4',alpha=0.5)


plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.3],[1.,1.],'C7--')
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$M^c_{200}/M_{200}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratioc_NMbin.pdf',bbox_inches='tight')

# -----------------------


plt.figure()

plt.scatter(lM200[mN1],pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23],pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],pcc[mN4M],facecolor='none',edgecolors='C8',label='$ 4 \leq N_{GAL} \leq 10$')

plt.errorbar(lM200[mN1],pcc[mN1],xerr=elM200[:,mN1],yerr=epcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23],pcc[mN23],xerr=elM200[:,mN23],yerr=epcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],pcc[mN4M],xerr=elM200[:,mN4M],yerr=epcc[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot(lM200c[mN1] ,pcc_c[mN1]  ,'C6^',label='$N_{GAL} = 1$')
plt.plot(lM200c[mN23],pcc_c[mN23],'C7^',label='$ 2 \leq N_{GAL} \leq 3$')
plt.plot(lM200c[mN4M],pcc_c[mN4M],'C8^',label='$ 4 \leq N_{GAL} \leq 10$')

plt.errorbar(lM200c[mN1] ,pcc_c[mN1] ,xerr=elM200c[:,mN1] ,yerr=epcc_c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200c[mN23],pcc_c[mN23],xerr=elM200c[:,mN23],yerr=epcc_c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200c[mN4M],pcc_c[mN4M],xerr=elM200c[:,mN4M],yerr=epcc_c[:,mN4M],fmt = 'none',ecolor='C8')


# plt.legend(frameon = False)
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/M200_pcc_NMbin.pdf',bbox_inches='tight')


# -----------------------
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
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratiopcc_NMbin.pdf',bbox_inches='tight')

# -----------------------
