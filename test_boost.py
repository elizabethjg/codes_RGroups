import sys
import numpy as np
from matplotlib import *
from astropy.io import fits
from matplotlib import rc
from astropy.cosmology import LambdaCDM
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 16})
plt.rc('font', family='serif', size='16.0')

cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

'''
test boost factor
'''

path = '../../'

#parameters
cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)


def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er

RIN = 300.
ROUT = 5000.
ndots = 15


bines = np.logspace(np.log10(RIN),np.log10(ROUT),num=ndots+1)
area = np.pi*np.diff(bines**2)

plt.figure()
for j in range(7):

    p = fits.open(path+'Rgroups/profiles_indcat/profile_Mbin'+str(j+1)+'.fits')[1].data

    ngal = p.NGAL_w

    d = ngal/area

    fcl = ((d - d[-1])*area)/ngal
    
    plt.plot(p.Rp,1./(1-fcl),label = 'Mbin'+str(j+1))
    
plt.legend()

plt.xlabel('R [Mpc]')
plt.ylabel('$1/(1-f_{cl})$')

plt.savefig(path+'Rgroups/plots_newanalysis/fcl.pdf',bbox_inches='tight')


out_T   = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_Mbin.out').T
out    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_nonboost.out').T

MH   = out[6]
MH_T = out_T[6]

lMH  = np.log10(out[6]*1.e14)
lMH_T  = np.log10(out_T[6]*1.e14)

M200  = out[11] 
M200_T  = out_T[11] 

eM200   = np.array([out[12],out[13]])
eM200_T   = np.array([out_T[12],out_T[13]])

lM200_T  = np.log10(out_T[11]*1.e14)
elM200_T = np.array([out_T[12]/(out_T[11]*np.log(10.)),out_T[13]/(out_T[11]*np.log(10.))])

lM200  = np.log10(out[11]*1.e14)
elM200   = np.array([out[12]/(out[11]*np.log(10.)),out[13]/(out[11]*np.log(10.))])

pcc   = out[14]
pcc_T = out_T[14]

epcc     = np.array([out[15],out[16]])
epcc_T     = np.array([out_T[15],out_T[16]])


# -----------------------
f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(lMH_T,lM200_T,facecolor='none',edgecolors='k'  ,label='Total sample')
ax[0].scatter(lMH,lM200,facecolor='none',edgecolors='C9'  ,label='non boost correc')

ax[0].errorbar(lMH_T,lM200_T,yerr=elM200_T,fmt = 'none',ecolor='k')
ax[0].errorbar(lMH,lM200,yerr=elM200,fmt = 'none',ecolor='C9')


ax[0].plot([12.3,15],[12.3,15],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_ylabel(r'$\log (M_{WL})$')
ax[0].axis([12.5,14.9,12.3,14.9])

ax[1].scatter(lMH_T  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
ax[1].scatter(lMH  ,M200/MH,facecolor='none',edgecolors='C9',label='Total sample')

ax[1].errorbar(lMH_T  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1].errorbar(lMH  ,M200/MH,yerr=eM200/MH,fmt = 'none',ecolor='C9')


ax[1].plot([12.3,15],[1.,1.],'C7--')
ax[1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')

plt.savefig(path+'Rgroups/plots_newanalysis/boost_correc2.pdf',bbox_inches='tight')

# -----------------------
f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(lM200_T,lM200,facecolor='none',edgecolors='k' )


ax[0].errorbar(lM200_T,lM200,xerr=elM200_T,yerr=elM200,fmt = 'none',ecolor='k')



ax[0].plot([12.0,15],[12.0,15],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_xlabel(r'$\log (M_{WL})$')
ax[0].set_ylabel(r'$\log (M^{wo}_{WL})$')
ax[0].axis([12.0,14.9,12.3,14.9])

ax[1].scatter(lM200_T  ,M200/M200_T,facecolor='none',edgecolors='k',label='Total sample')


ax[1].errorbar(lM200_T  ,M200/M200_T,yerr=eM200/M200_T,fmt = 'none',ecolor='k')



ax[1].plot([12.0,15],[1.,1.],'C7--')
ax[1].set_xlabel(r'$\log (M_{WL})$')
ax[1].set_ylabel(r'$M^{wo}_{WL}/M_{WL}$')

plt.savefig(path+'Rgroups/plots_newanalysis/boost_correc.pdf',bbox_inches='tight')
