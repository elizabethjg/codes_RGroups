import sys
import numpy as np
from astropy.io import fits
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 16})
plt.rc('font', family='serif', size='16.0')


out_T  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_Mbinb.out').T
out_H  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zHb.out').T
out_L  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zLb.out').T

out_Tc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_Mbinb_cM.out').T
out_Hc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zHb_cM.out').T
out_Lc  = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_zLb_cM.out').T


lMH_m  = np.log10(out_T[4]*1.e14)
lMH_M  = np.log10(out_T[5]*1.e14)

NL_T    = out_T[0]
NL_H    = out_H[0]
NL_L    = out_L[0]

NL_Tc   = out_Tc[0]
NL_Hc   = out_Hc[0]
NL_Lc   = out_Lc[0]

M200_T  = out_T[11]*10.
M200_H  = out_H[11]*10.
M200_L  = out_L[11]*10.

M200_Tc = out_Tc[11]*10.
M200_Hc = out_Hc[11]*10.
M200_Lc = out_Lc[11]*10.

eM200_T = np.array([out_T[12],out_T[13]])*10.
eM200_H = np.array([out_H[12],out_H[13]])*10.
eM200_L = np.array([out_L[12],out_L[13]])*10.

eM200_Tc = np.array([out_Tc[12],out_Tc[13]])*10.
eM200_Hc = np.array([out_Hc[12],out_Hc[13]])*10.
eM200_Lc = np.array([out_Lc[12],out_Lc[13]])*10.

pcc_T   = out_T[14]
pcc_H   = out_H[14]
pcc_L   = out_L[14]

pcc_Tc  = out_Tc[14]
pcc_Hc  = out_Hc[14]
pcc_Lc  = out_Lc[14]

epcc_T  = np.array([out_T[15],out_T[16]])
epcc_H  = np.array([out_H[15],out_H[16]])
epcc_L  = np.array([out_L[15],out_L[16]])

epcc_Tc = np.array([out_Tc[15],out_Tc[16]])
epcc_Hc = np.array([out_Hc[15],out_Hc[16]])
epcc_Lc = np.array([out_Lc[15],out_Lc[16]])

f=open('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/table_final.tex','w')

for j in range(len(NL_T)):
    
    f.write('& $['+str('%.1f' % (lMH_m[j]))+','+str('%.1f' % (lMH_M[j]))+')$ & ')

    f.write('$'+str(int(NL_T[j]))+'$ & ')
    if M200_T[j] < 3:
        f.write('$'+str('%.2f' % (M200_T[j]))+'_{-'+str('%.2f' % (eM200_T[0][j]))+'}^{+'+str('%.2f' % (eM200_T[1][j]))+'}$ & ')
    elif M200_T[j] < 10:
        f.write('$'+str('%.1f' % (M200_T[j]))+'_{-'+str('%.1f' % (eM200_T[0][j]))+'}^{+'+str('%.1f' % (eM200_T[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200_T[j]))+'_{-'+str(int(eM200_T[0][j]))+'}^{+'+str(int(eM200_T[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc_T[j]))+'_{-'+str('%.2f' % (epcc_T[0][j]))+'}^{+'+str('%.2f' % (epcc_T[1][j]))+'}$ & ')

    f.write('$'+str(int(NL_Tc[j]))+'$ & ')
    if M200_Tc[j] < 3:
        f.write('$'+str('%.2f' % (M200_Tc[j]))+'_{-'+str('%.2f' % (eM200_Tc[0][j]))+'}^{+'+str('%.2f' % (eM200_Tc[1][j]))+'}$ & ')
    elif M200_Tc[j] < 10:
        f.write('$'+str('%.1f' % (M200_Tc[j]))+'_{-'+str('%.1f' % (eM200_Tc[0][j]))+'}^{+'+str('%.1f' % (eM200_Tc[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200_Tc[j]))+'_{-'+str(int(eM200_Tc[0][j]))+'}^{+'+str(int(eM200_Tc[1][j]))+'}$ & ')

    f.write('$'+str('%.2f' % (pcc_Tc[j]))+'_{-'+str('%.2f' % (epcc_Tc[0][j]))+'}^{+'+str('%.2f' % (epcc_Tc[1][j]))+r'}$ \\ \vspace{0.15cm} ')
    f.write(' \n')

for j in range(len(NL_T)):
    
    f.write('& $['+str('%.1f' % (lMH_m[j]))+','+str('%.1f' % (lMH_M[j]))+')$ & ')
    
    f.write('$'+str(int(NL_L[j]))+'$ & ')
    if M200_L[j] < 3:
        f.write('$'+str('%.2f' % (M200_L[j]))+'_{-'+str('%.2f' % (eM200_L[0][j]))+'}^{+'+str('%.2f' % (eM200_L[1][j]))+'}$ & ')
    elif M200_L[j] < 10:
        f.write('$'+str('%.1f' % (M200_L[j]))+'_{-'+str('%.1f' % (eM200_L[0][j]))+'}^{+'+str('%.1f' % (eM200_L[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200_L[j]))+'_{-'+str(int(eM200_L[0][j]))+'}^{+'+str(int(eM200_L[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc_L[j]))+'_{-'+str('%.2f' % (epcc_L[0][j]))+'}^{+'+str('%.2f' % (epcc_L[1][j]))+'}$ & ')

    f.write('$'+str(int(NL_Lc[j]))+'$ & ')
    if M200_Lc[j] < 3:
        f.write('$'+str('%.2f' % (M200_Lc[j]))+'_{-'+str('%.2f' % (eM200_Lc[0][j]))+'}^{+'+str('%.2f' % (eM200_Lc[1][j]))+'}$ & ')
    elif M200_Lc[j] < 10:
        f.write('$'+str('%.1f' % (M200_Lc[j]))+'_{-'+str('%.1f' % (eM200_Lc[0][j]))+'}^{+'+str('%.1f' % (eM200_Lc[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200_Lc[j]))+'_{-'+str(int(eM200_Lc[0][j]))+'}^{+'+str(int(eM200_Lc[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc_Lc[j]))+'_{-'+str('%.2f' % (epcc_Lc[0][j]))+'}^{+'+str('%.2f' % (epcc_Lc[1][j]))+r'}$ \\ \vspace{0.15cm} ')

    f.write(' \n')

for j in range(len(NL_T)):
    
  
    f.write('& $['+str('%.1f' % (lMH_m[j]))+','+str('%.1f' % (lMH_M[j]))+')$ & ')

    f.write('$'+str(int(NL_H[j]))+'$ & ')
    if M200_H[j] < 3:
        f.write('$'+str('%.2f' % (M200_H[j]))+'_{-'+str('%.2f' % (eM200_H[0][j]))+'}^{+'+str('%.2f' % (eM200_H[1][j]))+'}$ & ')
    elif M200_H[j] < 10:
        f.write('$'+str('%.1f' % (M200_H[j]))+'_{-'+str('%.1f' % (eM200_H[0][j]))+'}^{+'+str('%.1f' % (eM200_H[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200_H[j]))+'_{-'+str(int(eM200_H[0][j]))+'}^{+'+str(int(eM200_H[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc_H[j]))+'_{-'+str('%.2f' % (epcc_H[0][j]))+'}^{+'+str('%.2f' % (epcc_H[1][j]))+r'}$ &')
    

    f.write('$'+str(int(NL_Hc[j]))+'$ & ')
    if M200_Hc[j] < 3:
        f.write('$'+str('%.2f' % (M200_Hc[j]))+'_{-'+str('%.2f' % (eM200_Hc[0][j]))+'}^{+'+str('%.2f' % (eM200_Hc[1][j]))+'}$ & ')
    elif M200_Hc[j] < 10:
        f.write('$'+str('%.1f' % (M200_Hc[j]))+'_{-'+str('%.1f' % (eM200_Hc[0][j]))+'}^{+'+str('%.1f' % (eM200_Hc[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200_Hc[j]))+'_{-'+str(int(eM200_Hc[0][j]))+'}^{+'+str(int(eM200_Hc[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc_Hc[j]))+'_{-'+str('%.2f' % (epcc_Hc[0][j]))+'}^{+'+str('%.2f' % (epcc_Hc[1][j]))+r'}$ \\ \vspace{0.15cm} ')


    f.write(' \n')
    

f.close()

out    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_Nbin.out').T
outc    = np.loadtxt('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/Lens_Nbin_cM.out').T



lMH_m  = np.log10(out[4]*1.e14)
lMH_M  = np.log10(out[5]*1.e14)

NL    = out[0]
NLc   = outc[0]

M200  = out[11]*10.
M200c = outc[11]*10.

eM200   = np.array([out[12],out[13]])*10.
eM200c  = np.array([outc[12],outc[13]])*10.

pcc   = out[14]
pcc_c = outc[14]

epcc   = np.array([outc[15],outc[16]])
epcc_c = np.array([outc[15],outc[16]])

f=open('/home/eli/Documentos/Astronomia/posdoc/Rgroups/profiles_newanalysis/table_final_Nbin.tex','w')

for j in range(len(NL)):
  
    f.write('& $['+str('%.1f' % (lMH_m[j]))+','+str('%.1f' % (lMH_M[j]))+')$ & ')

    f.write('$'+str(int(NL[j]))+'$ & ')
    if M200[j] < 3:
        f.write('$'+str('%.2f' % (M200[j]))+'_{-'+str('%.2f' % (eM200[0][j]))+'}^{+'+str('%.2f' % (eM200[1][j]))+'}$ & ')
    elif M200[j] < 10:
        f.write('$'+str('%.1f' % (M200[j]))+'_{-'+str('%.1f' % (eM200[0][j]))+'}^{+'+str('%.1f' % (eM200[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200[j]))+'_{-'+str(int(eM200[0][j]))+'}^{+'+str(int(eM200[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc[j]))+'_{-'+str('%.2f' % (epcc[0][j]))+'}^{+'+str('%.2f' % (epcc[1][j]))+r'}$ &')
    

    f.write('$'+str(int(NLc[j]))+'$ & ')
    if M200c[j] < 3:
        f.write('$'+str('%.2f' % (M200c[j]))+'_{-'+str('%.2f' % (eM200c[0][j]))+'}^{+'+str('%.2f' % (eM200c[1][j]))+'}$ & ')
    elif M200c[j] < 10:
        f.write('$'+str('%.1f' % (M200c[j]))+'_{-'+str('%.1f' % (eM200c[0][j]))+'}^{+'+str('%.1f' % (eM200c[1][j]))+'}$ & ')
    else:
        f.write('$'+str(int(M200c[j]))+'_{-'+str(int(eM200c[0][j]))+'}^{+'+str(int(eM200c[1][j]))+'}$ & ')
    f.write('$'+str('%.2f' % (pcc_c[j]))+'_{-'+str('%.2f' % (epcc_c[0][j]))+'}^{+'+str('%.2f' % (epcc_c[1][j]))+r'}$ \\ \vspace{0.15cm} ')
    
    f.write(' \n')

f.close()
