import numpy as np
# from pylab import *
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 12})
cosmo = LambdaCDM(H0=70., Om0=0.3, Ode0=0.7)

folder =  '/home/eli/Documentos/Astronomia/posdoc/Rgroups/haloellip/'
# folder = '/home/eli/Documentos/PostDoc/halo-elongation/redMapper/'

angles        = fits.open(folder+'angle_Rgroups_FINAL.fits')[1].data
coordinates   = fits.open(folder+'RGroup_projected_member_position.fits')[1].data
cat           = fits.open('/home/eli/Documentos/Astronomia/posdoc/Rgroups/catalogs/FINAL_Group.fits')[1].data


ind = np.argsort(angles.ID)
angles = angles[ind]

ID_c      = angles.ID
ID_c_tot  = cat.ID

mid = np.in1d(ID_c_tot,ID_c)
cat = cat[mid]

mid = np.in1d(coordinates.ID,ID_c)
coordinates = coordinates[mid]

ides  = coordinates.ID
ID,c = np.unique(ides,return_counts=True)

zc = cat.Z
D_ang    = np.array(cosmo.angular_diameter_distance(zc))
kpcscale = D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
KPCSCALE = np.repeat(kpcscale,c)
Rhalo    = np.repeat(cat.RADIUS_HALO,c)

t     = np.repeat(angles['theta'],c)
e     = np.repeat(angles['e'],c)

dx = (coordinates['X']*KPCSCALE*1.e-3)/Rhalo
dy = (coordinates['Y']*KPCSCALE*1.e-3)/Rhalo

mask = np.isfinite(dx)
dx = dx[mask]
dy = dy[mask]

x_t  = (dx*np.cos(t[mask]) + dy*np.sin(t[mask])) 
y_t = (-1.*dx*np.sin(t[mask]) + dy*np.cos(t[mask])) 



lgrid = 20
Max = 1800


xedges = np.linspace(-2.,2.,30)
yedges = np.linspace(-2.,2.,30)
xcenters = (xedges[:-1] + xedges[1:]) / 2.
ycenters = (yedges[:-1] + yedges[1:]) / 2.
X,Y = np.meshgrid(xcenters,ycenters)
f, ax = plt.subplots(1, 2, sharey=True,figsize=(8,4))


levels = np.linspace(10,200,5)
H, xedges, yedges = np.histogram2d(dx, dy, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[0].contour(X, Y, H.T,levels,cmap='plasma')


H, xedges, yedges = np.histogram2d(x_t,y_t, bins=(xedges, yedges))#,weights=np.log10(Lum_r))
ax[1].contour(X, Y, H.T,levels,cmap='plasma')
f.subplots_adjust(hspace=0,wspace=0)

# plt.setp(ax[0,1].get_xticklabels(), visible=False)


# ax[0,0].axis([-1.2,1.2,-1.2,1.2])

