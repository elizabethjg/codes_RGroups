import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/python_codes')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
import numpy as np
from pylab import *
from multipoles_shear import *
import emcee
import time
from multiprocessing import Pool
import argparse
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from profiles_fit import *
#parameters

cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)


parser = argparse.ArgumentParser()
parser.add_argument('-folder', action='store', dest='folder',default='./')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ncores', action='store', dest='ncores', default=15)
parser.add_argument('-Yanmiss', action='store', dest='yan', default='False')
args = parser.parse_args()

h_cosmo = 1.0

folder    = args.folder
file_name = args.file_name
ncores    = args.ncores
ncores    = int(ncores)

if 'True' in args.yan:
    ymiss      = True
    outfile = 'monopole_pcconly_ymiss_'+file_name[:-4]+'out'
    tau = 0.3
elif 'False' in args.yan:
    ymiss     = False
    outfile = 'monopole_pcconly_'+file_name[:-4]+'out'
    tau = 0.4



print 'fitting monopole misscentred'
print folder
print file_name

profile = fits.open(folder+file_name)
h       = profile[1].header
p       = profile[1].data
zmean   = h['Z_MEAN']    
Mhalo   = 10**h['lMASS_HALO_mean']
Rmean   = h['RADIUS_HALO_mean']
ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))/h_cosmo
soff = tau*Rmean



def log_likelihood(data_model, r, Gamma, e_Gamma):
    log_M200,pcc = data_model
    M200 = 10**log_M200
    multipoles = multipole_shear(r,M200=M200,
                                misscentred = True,s_off = soff,h = h_cosmo,
                                ellip=0,z=zmean,components = ['t'],
                                verbose=False,Yanmiss=ymiss) #,ncores=ncores)
    model = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc)
    sigma2 = e_Gamma**2
    return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
    

def log_probability(data_model, r, Gamma, e_Gamma):
    log_M200, pcc = data_model
    if 11. < log_M200 < 16.5 and 0.2 < pcc < 1.0:
        return log_likelihood(data_model, r, Gamma, e_Gamma)
    return -np.inf

# initializing

pos = np.array([np.random.uniform(11.5,16.0,10),
                np.random.uniform(0.3,1.0,10)]).T


nwalkers, ndim = pos.shape

#-------------------
# running emcee

pool = Pool(processes=(ncores))

maskr = (p.Rp < ROUT)

t1 = time.time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp[maskr],p.DSigma_T[maskr],p.error_DSigma_T[maskr]),
                                pool = pool)
sampler.run_mcmc(pos, 100, progress=True)
print '//////////////////////'
print '         TIME         '
print '----------------------'
print (time.time()-t1)/60.
pool.terminate()
#-------------------

#flat_samples = sampler.get_chain(discard=100, flat=True)

#p1 = np.percentile(flat_samples[:, 0], [16, 50, 84])

#print p1[1],np.diff(p1)

# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)

f1=open(folder+outfile,'w')
f1.write('# log(M200)  pcc  \n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f']*2)
f1.close()
