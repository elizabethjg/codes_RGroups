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

out_CS82 = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_indcat/Lens_CS82.out').T
out_CFHT = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_indcat/Lens_CFHT.out').T
out_KiDS = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_indcat/Lens_KiDS.out').T
out_RCSL = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_indcat/Lens_RCSL.out').T

lMH_CS82  = np.log10(out_CS82[6]*1.e14)
lMH_CFHT  = np.log10(out_CFHT[6]*1.e14)
lMH_KiDS  = np.log10(out_KiDS[6]*1.e14)
lMH_RCSL  = np.log10(out_RCSL[6]*1.e14)

M200_CS82  = out_CS82[11] 
M200_CFHT  = out_CFHT[11] 
M200_KiDS  = out_KiDS[11] 
M200_RCSL  = out_RCSL[11] 

eM200_CS82   = [out_CS82[12],out_CS82[13]]
eM200_CFHT   = [out_CFHT[12],out_CFHT[13]]
eM200_KiDS   = [out_KiDS[12],out_KiDS[13]]
eM200_RCSL   = [out_RCSL[12],out_RCSL[13]]

lM200_CS82  = np.log10(out_CS82[11]*1.e14)
lM200_CFHT  = np.log10(out_CFHT[11]*1.e14)
lM200_KiDS  = np.log10(out_KiDS[11]*1.e14)
lM200_RCSL  = np.log10(out_RCSL[11]*1.e14)

elM200_CS82   = [out_CS82[12]/(out_CS82[11]*np.log(10.)),out_CS82[13]/(out_CS82[11]*np.log(10.))]
elM200_CFHT   = [out_CFHT[12]/(out_CFHT[11]*np.log(10.)),out_CFHT[13]/(out_CFHT[11]*np.log(10.))]
elM200_KiDS   = [out_KiDS[12]/(out_KiDS[11]*np.log(10.)),out_KiDS[13]/(out_KiDS[11]*np.log(10.))]
elM200_RCSL   = [out_RCSL[12]/(out_RCSL[11]*np.log(10.)),out_RCSL[13]/(out_RCSL[11]*np.log(10.))]

pcc_CS82 = out_CS82[14]
pcc_CFHT = out_CFHT[14]
pcc_KiDS = out_KiDS[14]
pcc_RCSL = out_RCSL[14]

epcc_CS82 = np.array([out_CS82[15],out_CS82[16]])
epcc_CFHT = np.array([out_CFHT[15],out_CFHT[16]])
epcc_KiDS = np.array([out_KiDS[15],out_KiDS[16]])
epcc_RCSL = np.array([out_RCSL[15],out_RCSL[16]])

# -----------------------
plt.figure()

plt.scatter(lMH_CS82,lM200_CS82,facecolor='none',edgecolors='C0', label = 'CS82')
plt.scatter(lMH_CFHT,lM200_CFHT,facecolor='none',edgecolors='C1', label = 'CFHT')
plt.scatter(lMH_KiDS,lM200_KiDS,facecolor='none',edgecolors='C2', label = 'KiDS')
plt.scatter(lMH_RCSL,lM200_RCSL,facecolor='none',edgecolors='C3', label = 'RCSL')

plt.errorbar(lMH_CS82,lM200_CS82,yerr=elM200_CS82,fmt = 'none',ecolor='C0')
plt.errorbar(lMH_CFHT,lM200_CFHT,yerr=elM200_CFHT,fmt = 'none',ecolor='C1')
plt.errorbar(lMH_KiDS,lM200_KiDS,yerr=elM200_KiDS,fmt = 'none',ecolor='C2')
plt.errorbar(lMH_RCSL,lM200_RCSL,yerr=elM200_RCSL,fmt = 'none',ecolor='C3')

plt.plot([12.,15],[12.,15],'C7--')
plt.legend(frameon = False,fontsize=14)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{200})$')
plt.axis([12.,15,12.,15])
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_M200_Mbin_indcat.pdf',bbox_inches='tight')
# -----------------------

plt.figure()

plt.scatter(lM200_CS82,pcc_CS82,facecolor='none',edgecolors='C0', label = 'CS82')
plt.scatter(lM200_CFHT,pcc_CFHT,facecolor='none',edgecolors='C1', label = 'CFHT')
plt.scatter(lM200_KiDS,pcc_KiDS,facecolor='none',edgecolors='C2', label = 'KiDS')
plt.scatter(lM200_RCSL,pcc_RCSL,facecolor='none',edgecolors='C3', label = 'RCSL')

plt.errorbar(lM200_CS82,pcc_CS82,xerr=elM200_CS82,yerr=epcc_CS82,fmt = 'none',ecolor='C0')
plt.errorbar(lM200_CFHT,pcc_CFHT,xerr=elM200_CFHT,yerr=epcc_CFHT,fmt = 'none',ecolor='C1')
plt.errorbar(lM200_KiDS,pcc_KiDS,xerr=elM200_KiDS,yerr=epcc_KiDS,fmt = 'none',ecolor='C2')
plt.errorbar(lM200_RCSL,pcc_RCSL,xerr=elM200_RCSL,yerr=epcc_RCSL,fmt = 'none',ecolor='C3')


# plt.legend(frameon = False)
plt.xlabel('$\log (M_{200})$')
plt.ylabel('$p_{cc}$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/M200_pcc_Mbin_indcat.pdf',bbox_inches='tight')


