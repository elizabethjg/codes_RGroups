import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from astropy.wcs import WCS
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

ID_c  = groups.ID[mgal]

e = np.array([])
theta = np.array([])

X = np.array([])
Y = np.array([])

for j in range(len(ID_c)):
     print j
     mid  = ides == ID_c[j]
     N = np.append(N,mid.sum())
     ra   = RA[mid]
     dec  = DEC[mid]
     mcen = R_cen[mid] == 0
     
     ra0  = ra[mcen]
     dec0 = dec[mcen]
     
     wcs.wcs.crval = [ra0[0],dec0[0]]
     dx, dy = wcs.wcs_world2pix(ra, dec, 0)

     X = np.append(X,dx)
     Y = np.append(Y,dy)


     dx = dx[~mcen]
     dy = dy[~mcen]
     
     # all members
     ellip, ang = momentos(dx,dy,np.ones(len(dx)))
     e = np.append(e,ellip)
     theta = np.append(theta,ang)
     

X = np.array(X)
Y = np.array(Y)

theta2           = theta
theta2[theta<0.] = np.pi + theta[theta<0.]

# '''
tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='N', format='D', array=N),
        fits.Column(name='e', format='D', array=e),
        fits.Column(name='theta', format='D', array=theta)])
        
tbhdu.writeto(folder+'angle_Rgroups_FINAL.fits',overwrite=True)        


tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='X', format='D', array=X),
        fits.Column(name='Y', format='D', array=Y)])
        
tbhdu.writeto(folder+'RGroup_projected_member_position.fits',overwrite=True)        

# '''
