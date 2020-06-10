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

folder = '/home/eli/Documentos/Astronomia/posdoc/Rgroups/catalogs/'

Lcat = fits.open(folder+'FINAL_Group_Expanded.fits')
IDused = fits.open(folder+'ID_FINALcat_used.fits')[1].data.ID

mid = np.in1d(Lcat[1].data.ID,IDused)

lMH = np.log10(Lcat[1].data.MASS_HALO)
N = Lcat[1].data.N_GAL
c = Lcat[1].data.C_BG

mN1 = N==1
mN23 = (N>1)*(N<4)
mN4M = (N>3)
mc = c > 2.73

plt.figure()

plt.hist(lMH[mid*(lMH > 12.5)],np.linspace(12.5,15,70),histtype='step',density=True      ,label='Total sample',color='k')
plt.hist(lMH[mid*(lMH > 12.5)*mN1],np.linspace(12.5,15,70),histtype='step',density=True  ,label='$N_{GAL} = 1$',color='C6')
plt.hist(lMH[mid*(lMH > 12.5)*mN23],np.linspace(12.5,15,70),histtype='step',density=True ,label='$ 2 \leq N_{GAL} \leq 3$',color='C7')
plt.hist(lMH[mid*(lMH > 12.5)*mN4M],np.linspace(12.5,15,70),histtype='step',density=True,label='$ N_{GAL} \geq 4$',color='C8')
plt.legend(loc=1,fontsize=14,ncol=2,columnspacing=0.2,frameon=False)

plt.hist(lMH[mc*mid*(lMH > 12.5)],np.linspace(12.5,15,70),histtype='step',density=True,label='Total sample',color='k',ls='--')
plt.hist(lMH[mc*mid*(lMH > 12.5)*mN1],np.linspace(12.5,15,70),histtype='step',density=True,label='$N_{GAL} = 1$',color='C6',ls='--')
plt.hist(lMH[mc*mid*(lMH > 12.5)*mN23],np.linspace(12.5,15,70),histtype='step',density=True,label='$ 2 \leq N_{GAL} \leq 3$',color='C7',ls='--')
plt.hist(lMH[mc*mid*(lMH > 12.5)*mN4M],np.linspace(12.5,15,70),histtype='step',density=True,label='$ N_{GAL} \geq 4$',color='C8',ls='--')

plt.axis([12.5,15.,0.,2.3])
plt.xlabel('$\log{M_{AM}}$')
plt.ylabel('$n$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_dist_norm.pdf',bbox_inches='tight')

plt.figure()

plt.hist(lMH[mid*(lMH > 12.5)],np.linspace(12.5,15,70),histtype='step',density=False      ,label='Total sample',color='k')
plt.hist(lMH[mid*(lMH > 12.5)*mN1],np.linspace(12.5,15,70),histtype='step',density=False  ,label='$N_{GAL} = 1$',color='C6')
plt.hist(lMH[mid*(lMH > 12.5)*mN23],np.linspace(12.5,15,70),histtype='step',density=False ,label='$ 2 \leq N_{GAL} \leq 3$',color='C7')
plt.hist(lMH[mid*(lMH > 12.5)*mN4M],np.linspace(12.5,15,70),histtype='step',density=False,label='$ N_{GAL} \geq 4$',color='C8')
plt.legend(loc=1,fontsize=13,frameon=False)

plt.hist(lMH[mc*mid*(lMH > 12.5)],np.linspace(12.5,15,70),histtype='step',density=False,label='Total sample',color='k',ls='--')
plt.hist(lMH[mc*mid*(lMH > 12.5)*mN1],np.linspace(12.5,15,70),histtype='step',density=False,label='$N_{GAL} = 1$',color='C6',ls='--')
plt.hist(lMH[mc*mid*(lMH > 12.5)*mN23],np.linspace(12.5,15,70),histtype='step',density=False,label='$ 2 \leq N_{GAL} \leq 3$',color='C7',ls='--')
plt.hist(lMH[mc*mid*(lMH > 12.5)*mN4M],np.linspace(12.5,15,70),histtype='step',density=False,label='$ N_{GAL} \geq 4$',color='C8',ls='--')

plt.yscale('log')
# plt.axis([12.5,15.,0.,2.3])
plt.xlabel('$\log{M_{AM}}$')
plt.ylabel('$N$')
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/Mhalo_dist.pdf',bbox_inches='tight')



plt.figure()

plt.hist(c[mid*(lMH > 12.5)],np.linspace(1.5,4,70),histtype='step',density=True,label='Total sample',color='k')
plt.hist(c[mid*(lMH > 12.5)*mN1],np.linspace(1.5,4,70),histtype='step',density=True,label='$N_{GAL} = 1$',color='C6')
plt.hist(c[mid*(lMH > 12.5)*mN23],np.linspace(1.5,4,70),histtype='step',density=True,label='$ 2 \leq N_{GAL} \leq 3$',color='C7')
plt.hist(c[mid*(lMH > 12.5)*mN4M],np.linspace(1.5,4,70),histtype='step',density=True,label='$ N_{GAL} \geq 4$',color='C8')
plt.axvline(2.73,color='k',ls=':')

plt.axis([1.5,4.,0.,1.9])
plt.xlabel('$C$')
plt.ylabel('$n$')
# plt.legend(loc=2,fontsize=12,ncol=2,columnspacing=0.2,frameon=False)
plt.savefig('/home/eli/Documentos/Astronomia/posdoc/Rgroups/plots/con_dist.pdf',bbox_inches='tight')
