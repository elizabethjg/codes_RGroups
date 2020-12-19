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
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=100., Om0=0.3, Ode0=0.7)

folder = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/catalogs/'

Lcat = fits.open(folder+'FINAL_Group_Expanded.fits')
IDused = fits.open(folder+'ID_FINALcat_used.fits')[1].data.ID
Gcat = fits.open(folder+'FINAL_gals.fits')

mid = np.in1d(Lcat[1].data.ID,IDused)

z = Lcat[1].data.Z[mid]
ID = Lcat[1].data.ID[mid]


midg   = np.in1d(Gcat[1].data.group_id,IDused)
mag_r = Gcat[1].data.sdss_mag_r[midg]
IDg   = Gcat[1].data.group_id[midg]

L_r = np.array([])

for j in range(len(ID)):
    zc = z[j]
    mi = IDg == ID[j]
    D_lum    = np.array(cosmo.luminosity_distance(zc))*1.e6
    MAG_abs  = mag_r[mi] + 5.0 - 5.0*np.log10(D_lum)
    Lum_r    = 10.**(-0.4*MAG_abs)
    L_r = np.append(L_r,Lum_r.sum())

M_r = -2.5*np.log10(L_r)
lMH = np.log10(Lcat[1].data.MASS_HALO)[mid]
N = Lcat[1].data.N_GAL[mid]
c = Lcat[1].data.C_BG[mid]

Mmean = np.array([])
LH    = np.array([])
LL    = np.array([])
eLH   = np.array([])
eLL   = np.array([])

Mmeanc = np.array([])
LHc    = np.array([])
LLc    = np.array([])
eLHc   = np.array([])
eLLc   = np.array([])

Mmean1 = np.array([])
LH1    = np.array([])
LL1    = np.array([])
eLH1   = np.array([])
eLL1   = np.array([])

Mmean23 = np.array([])
LH23    = np.array([])
LL23    = np.array([])
eLH23   = np.array([])
eLL23   = np.array([])

Mmean410 = np.array([])
LH410    = np.array([])
LL410    = np.array([])
eLH410   = np.array([])
eLL410   = np.array([])

Mmean11M = np.array([])
LH11M    = np.array([])
LL11M    = np.array([])
eLH11M   = np.array([])
eLL11M   = np.array([])

mN1 = N==1
mN23 = (N>1)*(N<4)
mN410 = (N>3)*(N<11)
mN11M = N>10
mc    = c > 2.73

for j in range(10):
    mmin = 12.5+0.2*j
    print mmin
    mmass = (lMH >= mmin)*(lMH < (mmin + 0.2))

    mz    = z >= np.median(z[mmass]) 

    num = sum(mz*mmass)
    
    Mmean = np.append(Mmean,np.median(lMH[mz*mmass]))
    LH    = np.append(LH,np.median(L_r[mz*mmass]))
    LL    = np.append(LL,np.median(L_r[(~mz)*mmass]))    
    eLH    = np.append(eLH,np.std(L_r[mz*mmass])/np.sqrt(num))
    eLL    = np.append(eLL,np.std(L_r[(~mz)*mmass])/np.sqrt(num))    
    
    mask = mmass*mc
    mz    = z >= np.median(z[mask]) 
    num = sum(mask*mz)
    
    if num > 1:
        Mmeanc = np.append(Mmeanc,np.median(lMH[mask]))
        LHc    = np.append(LHc,np.median(L_r[mz*mask]))
        LLc    = np.append(LLc,np.median(L_r[(~mz)*mask]))    
        eLHc   = np.append(eLHc,np.std(L_r[mz*mask])/np.sqrt(num))
        eLLc   = np.append(eLLc,np.std(L_r[(~mz)*mask])/np.sqrt(num))    

    mask = mmass*mc*mN1
    mz    = z >= np.median(z[mask]) 
    num = sum(mask*mz)
    
    if num > 1:
        Mmean1 = np.append(Mmean1,np.median(lMH[mask]))
        LH1    = np.append(LH1,np.median(L_r[mz*mask]))
        LL1    = np.append(LL1,np.median(L_r[(~mz)*mask]))    
        eLH1   = np.append(eLH1,np.std(L_r[mz*mask])/np.sqrt(num))
        eLL1   = np.append(eLL1,np.std(L_r[(~mz)*mask])/np.sqrt(num))    

    mask = mmass*mc*mN23
    mz    = z >= np.median(z[mask]) 
    num = sum(mask*mz)
    
    if num > 1:
        Mmean23 = np.append(Mmean23,np.median(lMH[mask]))
        LH23    = np.append(LH23,np.median(L_r[mz*mask]))
        LL23    = np.append(LL23,np.median(L_r[(~mz)*mask]))    
        eLH23   = np.append(eLH23,np.std(L_r[mz*mask])/np.sqrt(num))
        eLL23   = np.append(eLL23,np.std(L_r[(~mz)*mask])/np.sqrt(num))    

    mask = mmass*mc*mN410
    mz    = z >= np.median(z[mask]) 
    num = sum(mask*mz)
    
    if num > 1:
        Mmean410 = np.append(Mmean410,np.median(lMH[mask]))
        LH410    = np.append(LH410,np.median(L_r[mz*mask]))
        LL410    = np.append(LL410,np.median(L_r[(~mz)*mask]))    
        eLH410   = np.append(eLH410,np.std(L_r[mz*mask])/np.sqrt(num))
        eLL410   = np.append(eLL410,np.std(L_r[(~mz)*mask])/np.sqrt(num))    

    mask = mmass*mc*mN11M
    mz    = z >= np.median(z[mask]) 
    num = sum(mask*mz)
    
    if num > 1:
        Mmean11M = np.append(Mmean11M,np.median(lMH[mask]))
        LH11M    = np.append(LH11M,np.median(L_r[mz*mask]))
        LL11M    = np.append(LL11M,np.median(L_r[(~mz)*mask]))    
        eLH11M   = np.append(eLH11M,np.std(L_r[mz*mask])/np.sqrt(num))
        eLL11M   = np.append(eLL11M,np.std(L_r[(~mz)*mask])/np.sqrt(num))    

    
ratio_LH  = LH/LL
eratio_LH = np.sqrt((eLH/LL)**2 + ((eLL*LH)/(LL**2))**2)

ratio_LHc  = LHc/LLc
eratio_LHc = np.sqrt((eLHc/LLc)**2 + ((eLLc*LHc)/(LLc**2))**2)

ratio_LH1  = LH1/LL1
eratio_LH1 = np.sqrt((eLH1/LL1)**2 + ((eLL1*LH1)/(LL1**2))**2)

ratio_LH23  = LH23/LL23
eratio_LH23 = np.sqrt((eLH23/LL23)**2 + ((eLL23*LH23)/(LL23**2))**2)

ratio_LH410  = LH410/LL410
eratio_LH410 = np.sqrt((eLH410/LL410)**2 + ((eLL410*LH410)/(LL410**2))**2)

ratio_LH11M  = LH11M/LL11M
eratio_LH11M = np.sqrt((eLH11M/LL11M)**2 + ((eLL11M*LH11M)/(LL11M**2))**2)


plt.figure()
plt.scatter(Mmean,ratio_LH,facecolor='none',edgecolors='k', label = 'Total sample')
plt.errorbar(Mmean,ratio_LH,yerr=eratio_LH,fmt = 'none',ecolor='k')

plt.scatter(Mmeanc,ratio_LHc,facecolor='none',edgecolors='C4', label = 'C-sample')
plt.errorbar(Mmeanc,ratio_LHc,yerr=eratio_LHc,fmt = 'none',ecolor='C4')
plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$L^H/L^L$')
plt.legend()
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratio_lum.pdf',bbox_inches='tight')


plt.figure()
plt.scatter(Mmean1,ratio_LH1,facecolor='none',edgecolors='C6', label = 'N = 1')
plt.errorbar(Mmean1,ratio_LH1,yerr=eratio_LH1,fmt = 'none',ecolor='C6')

plt.scatter(Mmean23,ratio_LH23,facecolor='none',edgecolors='C7', label = 'N = 2-3')
plt.errorbar(Mmean23,ratio_LH23,yerr=eratio_LH23,fmt = 'none',ecolor='C7')

plt.scatter(Mmean410,ratio_LH410,facecolor='none',edgecolors='C8', label = 'N = 4-10')
plt.errorbar(Mmean410,ratio_LH410,yerr=eratio_LH410,fmt = 'none',ecolor='C8')

plt.scatter(Mmean11M,ratio_LH11M,facecolor='none',edgecolors='C9', label = 'N $>$ 10')
plt.errorbar(Mmean11M,ratio_LH11M,yerr=eratio_LH11M,fmt = 'none',ecolor='C9')


plt.plot([12.3,14.6],[1.,1.],'C7--')
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$L^H/L^L$')
plt.legend(fontsize = 12)
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/ratio_lum_NMbin.pdf',bbox_inches='tight')


# plt.yscale('log')
