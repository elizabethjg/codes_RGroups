import sys
import time
import numpy as np
import pandas as pd
from lensing import LensCat, gentools
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)


with fits.open('FINAL_Group.fits') as f:
        ft = Table(f[1].data)
        cols = ['ID', 'RA', 'DEC', 'RA_BG', 'DEC_BG','Z', 'N_GAL', 'CONCENTRATION',
                'MASS_HALO', 'VDISP_HALO', 'RADIUS_HALO', 
                'MASS_DYN', 'VDISP_DYN', 'RADIUS_DYN']
        df = ft[cols].to_pandas()

        ra = df['RA']
        dec = df['DEC']
        in_cs82 = ((((ra<60.)&(dec<2.)&(dec>-2.)) | ((ra>300.)&(dec<2.)&(dec>-2.))))#&mask_z # CS82
        in_kids = (ra>125.)&(ra<240.)&(dec>-5.)&(dec<5.)#&mask_z  # KiDS North
        in_cfht = (ra>28.)&(ra<41.)&(dec>-13.)&(dec<-2.)  # CFHT W1
        in_cfht |= (ra>129.)&(ra<139.)&(dec>-8.)&(dec<1.)  # CFHT W2
        in_cfht |= (ra>207.)&(ra<222.)&(dec>50.)&(dec<59.)  # CFHT W3
        in_cfht |= (ra>329.)&(ra<337.)&(dec>-2.)&(dec<6.)  # CFHT W4
        in_rcsl = np.ones(len(df), dtype=bool)

        z_min, z_max = 0.05, 0.40
        folder = 'N_all/'
        lensname = 'RM_FINAL'
        mask_z = (df['Z']>=z_min)&(df['Z']<z_max)
        #print z_min, z_max, (in_cs82*mask_z).sum(), (in_kids*mask_z).sum(), (in_cfht*mask_z).sum(), (in_rcsl*mask_z).sum()
        in_cs82 &= mask_z
        in_kids &= mask_z
        in_cfht &= mask_z
        in_rcsl &= mask_z
        
        N_min, N_max = 1, 1
        mask_N = (df['N_GAL']>=N_min)#&(df['N_GAL']<=N_max)
        print N_min, N_max, (in_cs82*mask_N).sum(), (in_kids*mask_N).sum(), (in_cfht*mask_N).sum(), (in_rcsl*mask_N).sum()
        in_cs82 &= mask_N
        in_kids &= mask_N
        in_cfht &= mask_N
        in_rcsl &= mask_N


R_Mpc = 7.5
njobs = 1
#'''
t0 = time.time()
print '*** CS82 ***'
if in_cs82.sum()!=0:
        # Search neighbors around each lens
        R_deg = gentools.Mpc2deg(R_Mpc=R_Mpc, z=df['Z'], cosmo=cosmo)
        #columns = ['RAJ2000','DECJ2000','Z_B','e1','e2','m','weight','ODDS','fitclass','MASK', 'BPZ_LOW95']
        LensCat.CS82.load()
        mask = (LensCat.CS82.data['Z_B']>z_min)*(LensCat.CS82.data['Z_B']<1.3)
        #mask *= LensCat.CS82.data['ODDS']>0.4
        LensCat.CS82.data = LensCat.CS82.data[mask]

        cs82_L, cs82_S=LensCat.CS82.find_neighbors(centre=df[['RA','DEC']][in_cs82], 
		upper_radii=R_deg[in_cs82], append_data=df[in_cs82], compressed=True, njobs=njobs)

        print 'Tiempo', (time.time() - t0)/60.


        # Arreglar coordenadas..
        mALFA0 = cs82_L.data['RA']>275
        cs82_L.data['RA'][mALFA0] -= 360.
	mALFA0 = cs82_L.data['RA_BG']>275
	cs82_L.data['RA_BG'][mALFA0] -= 360.
        mALFA0 = cs82_S.data['RAJ2000']>275
        cs82_S.data['RAJ2000'][mALFA0] -= 360.

        cs82_L.write_to(folder+'gx_'+cs82_L.name+'_L_'+lensname+'.fits', overwrite=True)
	cs82_S.write_to(folder+'gx_'+cs82_S.name+'_S_'+lensname+'.fits', overwrite=True)

        LensCat.CS82.drop()
        del cs82_L, cs82_S
        print 'Tiempo', (time.time() - t0)/60.
#'''

#####################################################################################################
t0 = time.time()
print '*** KiDS ***'
if in_kids.sum()!=0:
        # Search neighbors around each lens
        R_deg = gentools.Mpc2deg(R_Mpc=R_Mpc, z=df['Z'], cosmo=cosmo)
        LensCat.KiDS.load()
        mask = (LensCat.KiDS.data['Z_B']>z_min)*(LensCat.KiDS.data['Z_B']<1.3)
        LensCat.KiDS.data = LensCat.KiDS.data[mask]
        kids_L, kids_S=LensCat.KiDS.find_neighbors(centre=df[['RA','DEC']][in_kids], 
                upper_radii=R_deg[in_kids], append_data=df[in_kids], compressed=True, njobs=njobs)

        print 'Tiempo', (time.time() - t0)/60.

        kids_L.write_to(folder+'gx_'+kids_L.name+'_L_'+lensname+'.fits', overwrite=True)
        kids_S.write_to(folder+'gx_'+kids_S.name+'_S_'+lensname+'.fits', overwrite=True)
        
        LensCat.KiDS.drop()
        del kids_L, kids_S
        print 'Tiempo', (time.time() - t0)/60.
#'''
#####################################################################################################
t0 = time.time()
print '*** CFHT ***'
if in_cfht.sum()!=0:
        # Search neighbors around each lens
        R_deg = gentools.Mpc2deg(R_Mpc=R_Mpc, z=df['Z'], cosmo=cosmo)
        LensCat.CFHT.load()
        mask = (LensCat.CFHT.data['Z_B']>z_min)*(LensCat.CFHT.data['Z_B']<1.3)
        LensCat.CFHT.data = LensCat.CFHT.data[mask]
        cfht_L, cfht_S=LensCat.CFHT.find_neighbors(centre=df[['RA','DEC']][in_cfht], 
                upper_radii=R_deg[in_cfht], append_data=df[in_cfht], compressed=True, njobs=njobs)

        print 'Tiempo', (time.time() - t0)/60.

        cfht_L.write_to(folder+'gx_'+cfht_L.name+'_L_'+lensname+'.fits', overwrite=True)
        cfht_S.write_to(folder+'gx_'+cfht_S.name+'_S_'+lensname+'.fits', overwrite=True)

        LensCat.CFHT.drop()
        del cfht_L, cfht_S
        print 'Tiempo', (time.time() - t0)/60.

#####################################################################################################
t0 = time.time()
print '*** RCSL ***'
if in_rcsl.sum()!=0:
        # Search neighbors around each lens
        R_deg = gentools.Mpc2deg(R_Mpc=R_Mpc, z=df['Z'], cosmo=cosmo)
        LensCat.RCSL.load()
        mask = (LensCat.RCSL.data['Z_B']>z_min)*(LensCat.RCSL.data['Z_B']<1.3)
        LensCat.RCSL.data = LensCat.RCSL.data[mask]
        rcsl_L, rcsl_S=LensCat.RCSL.find_neighbors(centre=df[['RA','DEC']][in_rcsl], 
                upper_radii=R_deg[in_rcsl], append_data=df[in_rcsl], compressed=True, njobs=njobs)

        print 'Tiempo', (time.time() - t0)/60.

        rcsl_L.write_to(folder+'gx_'+rcsl_L.name+'_L_'+lensname+'.fits', overwrite=True)
        rcsl_S.write_to(folder+'gx_'+rcsl_S.name+'_S_'+lensname+'.fits', overwrite=True)

        LensCat.RCSL.drop()
        del rcsl_L, rcsl_S
        print 'Tiempo', (time.time() - t0)/60.
