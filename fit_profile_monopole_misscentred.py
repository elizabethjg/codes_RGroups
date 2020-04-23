import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
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

parser = argparse.ArgumentParser()
parser.add_argument('-folder', action='store', dest='folder',default='./')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ncores', action='store', dest='ncores', default=4)
args = parser.parse_args()

folder    = args.folder
file_name = args.file_name
ncores    = args.ncores
ncores    = int(ncores)

print 'fitting monopole misscentred'
print folder
print file_name

profile = fits.open(folder+file_name)
h       = profile[1].header
p       = profile[1].data
zmean   = h['Z_MEAN']    
Mhalo   = 10**h['lMASS_HALO_mean']
Rmean   = h['RADIUS_HALO_mean']
ROUT = (2.5*(2.*(Mhalo/2.21e14)**0.75)**(1./3.))/0.7

def log_likelihood(data_model, r, Gamma, e_Gamma):
    log_M200,pcc,tau = data_model
    M200 = 10**log_M200
    s_off = tau*Rmean
    multipoles = multipole_shear_parallel(r,M200=M200,
                                misscentred = True,s_off = s_off,
                                ellip=0,z=zmean,components = ['t'],
                                verbose=False,ncores=ncores)
    model = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc)
    sigma2 = e_Gamma**2
    return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
    

def log_probability(data_model, r, Gamma, e_Gamma):
    log_M200, pcc, tau = data_model
    if 11. < log_M200 < 15.5 and 0.3 < pcc < 1.0 and 0.1 < tau < 0.6:
        return log_likelihood(data_model, r, Gamma, e_Gamma)
    return -np.inf

# initializing

pos = np.array([np.random.uniform(11.5,15.0,10),
                np.random.normal(0.8,0.3,10),
                np.random.normal(0.1,0.6,10)]).T

pccdist = pos[:,1]                
pos[pccdist > 1.,1] = 1.

nwalkers, ndim = pos.shape

#-------------------
# running emcee

#pool = Pool(processes=(ncores))

maskr = (p.Rp < ROUT)

t1 = time.time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp[maskr],p.DSigma_T[maskr],p.error_DSigma_T[maskr]))
sampler.run_mcmc(pos, 500, progress=True)
print '//////////////////////'
print '         TIME         '
print '----------------------'
print (time.time()-t1)/60.
#pool.terminate()
#-------------------

#flat_samples = sampler.get_chain(discard=100, flat=True)

#p1 = np.percentile(flat_samples[:, 0], [16, 50, 84])

#print p1[1],np.diff(p1)

# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)

f1=open(folder+'monopole_misscentred_'+file_name,'w')
f1.write('# log(M200)  pcc  tau\n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f']*3)
f1.close()


mcmc1 = (mcmc_out.T)[:,500:]
M200,pcc = np.median(mcmc1,axis=1)

f=open(folder+file_name,'a')
f.write('#'+str('%.2f' % (10**(M200)/1.e14))+'   \n')
f.write('#'+str('%.2f' % (pcc))+'   \n')
f.close()
