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

def log_likelihood(data_model, r, Gamma, e_Gamma):
    log_M200,pcc = data_model
    M200 = 10**log_M200
    multipoles = multipole_shear_parallel(r,M200=M200,misscentred = True,
                                ellip=0,z=zmean,components = ['t'],
                                verbose=False,ncores=ncores)
    model = model_Gamma(multipoles,'t', misscentred = True, pcc = pcc)
    sigma2 = e_Gamma**2
    return -0.5 * np.sum((Gamma - model)**2 / sigma2 + np.log(2.*np.pi*sigma2))
    

def log_probability(data_model, r, Gamma, e_Gamma):
    log_M200, pcc = data_model
    if 12.5 < log_M200 < 15.5 and 0.3 < pcc < 1.0:
        return log_likelihood(data_model, r, Gamma, e_Gamma)
    return -np.inf

# initializing

pos = np.array([np.random.uniform(12.5,15.5,10),
                np.random.normal(0.8,0.3,10)]).T

pccdist = pos[:,1]                
pos[pccdist > 1.,1] = 1.

nwalkers, ndim = pos.shape

#-------------------
# running emcee

#pool = Pool(processes=(ncores))


t1 = time.time()
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp,p.DSigma_T,error_DSigma_T))
sampler.run_mcmc(pos, 200, progress=True)
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
f1.write('# log(M200)  pcc \n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f']*2)
f1.close()


mcmc1 = (mcmc_out.T)[:,500:]
M200,pcc = np.median(mcmc1,axis=1)

f=open(folder+file_name,'a')
f.write('#'+str('%.2f' % (10**(M200)/1.e14))+'   \n')
f.write('#'+str('%.2f' % (pcc))+'   \n')
f.close()