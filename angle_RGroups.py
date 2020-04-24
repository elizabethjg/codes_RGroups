import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from astropy.wcs import WCS
import matplotlib.pyplot as plt
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

def momentos(dx,dy,w):
     
     Q11  = np.sum((dx**2)*w)/np.sum(w)
     Q22  = np.sum((dy**2)*w)/np.sum(w)
     Q12  = np.sum((dx*dy)*w)/np.sum(w)
     E1 = (Q11-Q22)/(Q11+Q22)
     E2 = (2.*Q12)/(Q11+Q22)
     e = np.sqrt(E1**2 + E2**2)
     theta = np.arctan2(E2,E1)/2.
     return e,theta


wcs = WCS(naxis=2)
wcs.wcs.crpix = [0., 0.]
wcs.wcs.cdelt = [1./3600., 1./3600.]
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]    


groups  = fits.open('/mnt/clemente/lensing/RodriguezGroups/N_all/gx_L_RM_FINAL.fits')[1].data
members = fits.open('/mnt/clemente/lensing/RodriguezGroups/FINAL_gals.fits')[1].data

ides = members.group_id

RA0  = groups.RA_BG 
RA0[RA0 > 275] = RA0[RA0>275] - 360.
DEC0 = groups.DEC_BG 

RA    = members.RA
RA[RA > 275] = RA[RA>275] - 360.
DEC   = members.DEC

N_GAL = groups.N_GAL
mgal  = N_GAL > 2

ID    = groups.ID
RA0   = RA0
DEC0  = DEC0

ID_c = np.arange(len(ID))[mgal]

e     = np.zeros(len(ID))
theta = np.zeros(len(ID))

X = np.zeros(len(RA))
Y = np.zeros(len(RA))

for j in ID_c:
     print j, N_GAL[j]
     mid  = ides == ID[j]

     ra   = RA[mid]
     dec  = DEC[mid]
     
     ra0  = RA0[j]
     dec0 = DEC0[j]
     
     wcs.wcs.crval = [ra0,dec0]
     dx, dy = wcs.wcs_world2pix(ra, dec, 0)

     X[mid] = dx
     Y[mid] = dy
     
     # all members
     e[j], theta[j] = momentos(dx,dy,np.ones(len(dx)))
     
     l = 0.5*(max(dx)-min(dx))
     plt.figure()
     plt.plot(dx,dy,'k.')
     plt.plot(0,0,'ro')
     plt.plot([0.,l*np.cos(theta[j])],[0,l*np.sin(theta[j])],label = np.str(np.rad2deg(theta[j])))
     plt.legend()
     plt.savefig('/mnt/clemente/lensing/RodriguezGroups/'+str(j)+'.png',format='png',bbox_inches='tight')


theta2           = theta
theta2[theta<0.] = np.pi + theta[theta<0.]

# '''
tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='ID', format='D', array=ID),
        fits.Column(name='e', format='D', array=e),
        fits.Column(name='theta', format='D', array=theta2)])
        
tbhdu.writeto('/mnt/clemente/lensing/RodriguezGroups/angle_Rgroups_FINAL.fits',overwrite=True)        


tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='X', format='D', array=X),
        fits.Column(name='Y', format='D', array=Y)])
        
tbhdu.writeto('/mnt/clemente/lensing/RodriguezGroups/RGroup_projected_member_position.fits',overwrite=True)        

# '''
