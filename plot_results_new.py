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
Este es el codigo con el que se generaron los plots finales del 
draft
'''

#parameters
cvel = 299792458;   # Speed of light (m.s-1)
G    = 6.670e-11;   # Gravitational constant (m3.kg-1.s-2)
pc   = 3.085678e16; # 1 pc (m)
Msun = 1.989e30 # Solar mass (kg)


def ratio(x,ex,y,ey):
    r  = x/y
    er = np.sqrt((ex/y)**2 + ((x*ey)/(y**2))**2)
    return r, er

path = '../../'

pcc_s   = (np.loadtxt(path+'Rgroups/results_simu/pcc.txt').T)
masa_s   = (np.loadtxt(path+'Rgroups/results_simu/puntos-masamockmasaab_mass.txt').T)
z_s   = (np.loadtxt(path+'Rgroups/results_simu/puntos-masamockmasaab_z.txt').T)

out_T   = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_Mbin.out').T
out    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_Nbin.out').T

# lMdyn,slMdyn, sdisp     = np.loadtxt(path+'Rgroups/profiles_new/median_dyn_cM').T
lMdyn,slMdyn, sdisp     = np.loadtxt(path+'Rgroups/profiles_new/median_dyn').T

lMH_simu2 = masa_s[3]
ratio_simu = masa_s[4]

lMH_simu = pcc_s[2]
pcc_simu = 1-(pcc_s[-1]/100.)

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


#------------ CON SAMPLES ------------

outc    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_Nbin_cM.out').T
outc_T  = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_Mbin_cM.out').T

lMHc  = np.log10(outc[6]*1.e14)
lMHc_T  = np.log10(outc_T[6]*1.e14)

MHc   = outc[6]
MHc_T = outc_T[6]


M200c  = outc[11] 
M200c_T  = outc_T[11] 

eM200c   = np.array([outc[12],outc[13]])
eM200c_T  = np.array([outc_T[12],outc_T[13]])

ratio_c, eratio_c = ratio(M200c,eM200c,M200,eM200)
ratio_cT, eratio_cT = ratio(M200c_T,eM200c_T,M200_T,eM200_T)


lM200c_T  = np.log10(outc_T[11]*1.e14)


lM200c  = np.log10(outc[11]*1.e14)


elM200c   = np.array([outc[12]/(outc[11]*np.log(10.)),outc[13]/(outc[11]*np.log(10.))])
elM200c_T = np.array([outc_T[12]/(outc_T[11]*np.log(10.)),outc_T[13]/(outc_T[11]*np.log(10.))])

pcc_cT = outc_T[14]
pcc_c   = outc[14]

epcc_cT     = np.array([outc_T[15],outc_T[16]])
epcc_c     = np.array([outc[15],outc[16]])

ratio_pcc,eratio_pcc     = ratio(pcc_c,epcc_c,pcc,epcc)
ratio_pcc_T,eratio_pcc_T = ratio(pcc_cT,epcc_cT,pcc_T,epcc_T)

ratio_pcc  = np.append(ratio_pcc,ratio_pcc_T)
eratio_pcc = np.append(eratio_pcc,eratio_pcc_T)

#------------ IND CATS ------------

out_cs82    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_CS82.out').T
out_cfht    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_CFHT.out').T
out_rcsl    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_RCSL.out').T
out_kids    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_KiDS.out').T


lMH_cs82  = np.log10(out_cs82[6]*1.e14)
lMH_cfht  = np.log10(out_cfht[6]*1.e14)
lMH_rcsl  = np.log10(out_rcsl[6]*1.e14)
lMH_kids  = np.log10(out_kids[6]*1.e14)


MH_cs82   = out_cs82[6]
MH_cfht   = out_cfht[6]
MH_rcsl   = out_rcsl[6]
MH_kids   = out_kids[6]

M200_cs82 = out_cs82[11] 
M200_cfht = out_cfht[11] 
M200_rcsl = out_rcsl[11] 
M200_kids = out_kids[11] 

eM200_cs82 = np.array([out_cs82[12],out_cs82[13]])
eM200_cfht = np.array([out_cfht[12],out_cfht[13]])
eM200_rcsl = np.array([out_rcsl[12],out_rcsl[13]])
eM200_kids = np.array([out_kids[12],out_kids[13]])

lM200_cs82 = np.log10(out_cs82[11]*1.e14)
lM200_cfht = np.log10(out_cfht[11]*1.e14)
lM200_kids = np.log10(out_kids[11]*1.e14)
lM200_rcsl = np.log10(out_rcsl[11]*1.e14)

elM200_cs82 = np.array([out_cs82[12]/(out_cs82[11]*np.log(10.)),out_cs82[13]/(out_cs82[11]*np.log(10.))])
elM200_cfht = np.array([out_cfht[12]/(out_cfht[11]*np.log(10.)),out_cfht[13]/(out_cfht[11]*np.log(10.))])
elM200_kids = np.array([out_kids[12]/(out_kids[11]*np.log(10.)),out_kids[13]/(out_kids[11]*np.log(10.))])
elM200_rcsl = np.array([out_rcsl[12]/(out_rcsl[11]*np.log(10.)),out_rcsl[13]/(out_rcsl[11]*np.log(10.))])

pcc_cs82 = out_cs82[14]
pcc_cfht = out_cfht[14]
pcc_kids = out_kids[14]
pcc_rcsl = out_rcsl[14]

epcc_cs82 = np.array([out_cs82[15],out_cs82[16]])
epcc_cfht = np.array([out_cfht[15],out_cfht[16]])
epcc_kids = np.array([out_kids[15],out_kids[16]])
epcc_rcsl = np.array([out_rcsl[15],out_rcsl[16]])


#------------ TEST CATS ------------

out_o4    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_odds4.out').T
out_o6    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_odds6.out').T
out_z2    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_zcut2.out').T
out_z3    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_zcut3.out').T


lMH_o4  = np.log10(out_o4[6]*1.e14)
lMH_o6  = np.log10(out_o6[6]*1.e14)
lMH_z2  = np.log10(out_z2[6]*1.e14)
lMH_z3  = np.log10(out_z3[6]*1.e14)


MH_o4   = out_o4[6]
MH_o6   = out_o6[6]
MH_z2   = out_z2[6]
MH_z3   = out_z3[6]

M200_o4 = out_o4[11] 
M200_o6 = out_o6[11] 
M200_z2 = out_z2[11] 
M200_z3 = out_z3[11] 

eM200_o4 = np.array([out_o4[12],out_o4[13]])
eM200_o6 = np.array([out_o6[12],out_o6[13]])
eM200_z2 = np.array([out_z2[12],out_z2[13]])
eM200_z3 = np.array([out_z3[12],out_z3[13]])

lM200_o4 = np.log10(out_o4[11]*1.e14)
lM200_o6 = np.log10(out_o6[11]*1.e14)
lM200_z2 = np.log10(out_z2[11]*1.e14)
lM200_z3 = np.log10(out_z3[11]*1.e14)

elM200_o4 = np.array([out_o4[12]/(out_o4[11]*np.log(10.)),out_o4[13]/(out_o4[11]*np.log(10.))])
elM200_o6 = np.array([out_o6[12]/(out_o6[11]*np.log(10.)),out_o6[13]/(out_o6[11]*np.log(10.))])
elM200_z2 = np.array([out_z2[12]/(out_z2[11]*np.log(10.)),out_z2[13]/(out_z2[11]*np.log(10.))])
elM200_z3 = np.array([out_z3[12]/(out_z3[11]*np.log(10.)),out_z3[13]/(out_z3[11]*np.log(10.))])

pcc_o4 = out_o4[14]
pcc_o6 = out_o6[14]
pcc_z2 = out_z2[14]
pcc_z3 = out_z3[14]

epcc_o4 = np.array([out_o4[15],out_o4[16]])
epcc_o6 = np.array([out_o6[15],out_o6[16]])
epcc_z2 = np.array([out_z2[15],out_z2[16]])
epcc_z3 = np.array([out_z3[15],out_z3[16]])

############ Z samples

out_H    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_zH.out').T
out_L    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_zL.out').T

lMH_H  = np.log10(out_H[6]*1.e14)
lMH_L  = np.log10(out_L[6]*1.e14)

MH_H   = out_H[6]
MH_L   = out_L[6]

M200_H  = out_H[11] 
M200_L  = out_L[11] 

eM200_H   = np.array([out_H[12],out_H[13]])
eM200_L   = np.array([out_L[12],out_L[13]])

lM200_H  = np.log10(out_H[11]*1.e14)
lM200_L  = np.log10(out_L[11]*1.e14)

elM200_H = np.array([out_H[12]/(out_H[11]*np.log(10.)),out_H[13]/(out_H[11]*np.log(10.))])
elM200_L = np.array([out_L[12]/(out_L[11]*np.log(10.)),out_L[13]/(out_L[11]*np.log(10.))])

pcc_H = out_H[14]
pcc_L = out_L[14]

epcc_H = np.array([out_H[15],out_H[16]])
epcc_L = np.array([out_L[15],out_L[16]])


out_HcM    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_zH_cM.out').T
out_LcM    = np.loadtxt(path+'Rgroups/profiles_indcat/Lens_zL_cM.out').T

lMH_Hc  = np.log10(out_HcM[6]*1.e14)
lMH_Lc  = np.log10(out_LcM[6]*1.e14)

lM200_HcM  = np.log10(out_HcM[11]*1.e14)
lM200_LcM  = np.log10(out_LcM[11]*1.e14)

elM200_HcM = np.array([out_HcM[12]/(out_HcM[11]*np.log(10.)),out_HcM[13]/(out_HcM[11]*np.log(10.))])
elM200_LcM = np.array([out_LcM[12]/(out_LcM[11]*np.log(10.)),out_LcM[13]/(out_LcM[11]*np.log(10.))])


MH_HcM   = out_HcM[6]
MH_LcM   = out_LcM[6]

M200_HcM  = out_HcM[11] 
M200_LcM  = out_LcM[11] 

eM200_HcM   = np.array([out_HcM[12],out_HcM[13]])
eM200_LcM   = np.array([out_LcM[12],out_LcM[13]])

pcc_Hc = out_HcM[14]
pcc_Lc = out_LcM[14]

epcc_Hc = np.array([out_HcM[15],out_HcM[16]])
epcc_Lc = np.array([out_LcM[15],out_LcM[16]])

ratio_cH, eratio_cH = ratio(M200_HcM,eM200_HcM,M200_H,eM200_H)
ratio_cL, eratio_cL = ratio(M200_LcM,eM200_LcM,M200_L,eM200_L)



#---------------------------------------------------
lMH2 = np.append(lMH,lMH_T)

N = out[3]
mN1 = N==1
mN23 = (N>1)*(N<4)
mN4M = (N>3)

H        = cosmo.H(out[9][mN4M]).value/(1.0e3*pc) #H at z_pair s-1
M200_SIS  = ((2.*(sdisp*1.e3)**3)/((50**0.5)*G*H))/(Msun)

lMdyn2  = np.log10(M200_SIS)

# sdisp = out[-1][mN4M]

# -----------------------
# INDCATS

# f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f, ax2 = plt.subplots(figsize=(6.2,5))
ax = [ax2]

f.subplots_adjust(hspace=0,wspace=0)


ax[0].scatter(lM200_T[:-1],lM200_cs82,facecolor='none',edgecolors='C0'  ,label='CS82')
ax[0].scatter(lM200_T,lM200_cfht,facecolor='none',edgecolors='C1'  ,label='CFHT')
ax[0].scatter(lM200_T,lM200_rcsl,facecolor='none',edgecolors='C2'  ,label='RCSL')
ax[0].scatter(lM200_T,lM200_kids,facecolor='none',edgecolors='C3'  ,label='KiDS')


ax[0].errorbar(lM200_T[:-1],lM200_cs82,xerr=elM200_T[:,:-1],yerr=elM200_cs82,facecolor='none',ecolor='C0'  ,fmt = 'none')
ax[0].errorbar(lM200_T,lM200_cfht,xerr=elM200_T,yerr=elM200_cfht,facecolor='none',ecolor='C1'  ,fmt = 'none')
ax[0].errorbar(lM200_T,lM200_rcsl,xerr=elM200_T,yerr=elM200_rcsl,facecolor='none',ecolor='C2'  ,fmt = 'none')
ax[0].errorbar(lM200_T,lM200_kids,xerr=elM200_T,yerr=elM200_kids,facecolor='none',ecolor='C3'  ,fmt = 'none')


ax[0].plot([12.0,16],[12.0,16],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_ylabel('$\log (M^{ind}_{WL})$')
ax[0].axis([12.3,14.9,12.0,16])

# ax[1].scatter(lM200_T[:-1]  ,(M200_cs82 - M200_T[:-1])/np.mean(eM200_cs82,axis=0),facecolor='none',edgecolors='C0')
# ax[1].scatter(lM200_T  ,(M200_cfht - M200_T)/np.mean(eM200_cfht,axis=0),facecolor='none',edgecolors='C1')
# ax[1].scatter(lM200_T  ,(M200_rcsl - M200_T)/np.mean(eM200_rcsl,axis=0),facecolor='none',edgecolors='C2')
# ax[1].scatter(lM200_T  ,(M200_kids - M200_T)/np.mean(eM200_kids,axis=0),facecolor='none',edgecolors='C3')




# ax[1].plot([12.3,15],[0.,0.],'C7--')
# ax[1].axis([12.3,14.9,-2.5,2.5])

ax[0].set_xlabel(r'$\log (\langle M_{WL} \rangle)$')
# ax[1].set_ylabel(r'$(M^{ind}_{WL} - M_{WL})/ \epsilon(M^{ind}_{WL})$')

plt.savefig(path+'Rgroups/plots_newanalysis/M200_indcomp.pdf',bbox_inches='tight')

# -----------------------
# -----------------------
# INDCATS

f, ax = plt.subplots(figsize=(6.2,5))


ax.scatter(pcc_T[:-1],pcc_cs82,facecolor='none',edgecolors='C0'  ,label='CS82')
ax.scatter(pcc_T,pcc_cfht,facecolor='none',edgecolors='C1'  ,label='CFHT')
ax.scatter(pcc_T,pcc_rcsl,facecolor='none',edgecolors='C2'  ,label='RCSL')
ax.scatter(pcc_T,pcc_kids,facecolor='none',edgecolors='C3'  ,label='KiDS')


ax.errorbar(pcc_T[:-1],pcc_cs82,xerr=epcc_T[:,:-1],yerr=epcc_cs82,facecolor='none',ecolor='C0'  ,fmt = 'none')
ax.errorbar(pcc_T,pcc_cfht,xerr=epcc_T,yerr=epcc_cfht,facecolor='none',ecolor='C1'  ,fmt = 'none')
ax.errorbar(pcc_T,pcc_rcsl,xerr=epcc_T,yerr=epcc_rcsl,facecolor='none',ecolor='C2'  ,fmt = 'none')
ax.errorbar(pcc_T,pcc_kids,xerr=epcc_T,yerr=epcc_kids,facecolor='none',ecolor='C3'  ,fmt = 'none')


ax.plot([0.2,1.0],[0.2,1.0],'C7--')
ax.legend(frameon = False,fontsize = 13)

ax.set_ylabel('$p^{ind}_{cc}$')
ax.set_xlabel('$p_{cc}$')

plt.savefig(path+'Rgroups/plots_newanalysis/pcc_indcomp.pdf',bbox_inches='tight')

# -----------------------

# -----------------------
# TESTCATS

# f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f, ax2 = plt.subplots(figsize=(6.2,5))
ax = [ax2]

f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(lM200_T,lM200_o4,facecolor='none',edgecolors='C0'  ,label='ODDS $> 0.4$')
ax[0].scatter(lM200_T,lM200_o6,facecolor='none',edgecolors='C1'  ,label='ODDS $> 0.6$')
ax[0].scatter(lM200_T,lM200_z2,facecolor='none',edgecolors='C2'  ,label='$z_{cut} = 0.2$')
ax[0].scatter(lM200_T,lM200_z3,facecolor='none',edgecolors='C3'  ,label='$z_{cut} = 0.3$')


ax[0].errorbar(lM200_T,lM200_o4,xerr=elM200_T,yerr=elM200_o4,facecolor='none',ecolor='C0'  ,fmt = 'none')
ax[0].errorbar(lM200_T,lM200_o6,xerr=elM200_T,yerr=elM200_o6,facecolor='none',ecolor='C1'  ,fmt = 'none')
ax[0].errorbar(lM200_T,lM200_z2,xerr=elM200_T,yerr=elM200_z2,facecolor='none',ecolor='C2'  ,fmt = 'none')
ax[0].errorbar(lM200_T,lM200_z3,xerr=elM200_T,yerr=elM200_z3,facecolor='none',ecolor='C3'  ,fmt = 'none')


ax[0].plot([12.0,16],[12.0,16],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_ylabel('$\log (M^{ind}_{WL})$')
ax[0].axis([12.3,14.9,12.0,16])

# ax[1].scatter(lM200_T  ,(M200_o4 - M200_T)/np.mean(eM200_o4,axis=0),facecolor='none',edgecolors='C0')
# ax[1].scatter(lM200_T  ,(M200_o6 - M200_T)/np.mean(eM200_o6,axis=0),facecolor='none',edgecolors='C1')
# ax[1].scatter(lM200_T  ,(M200_z2 - M200_T)/np.mean(eM200_z2,axis=0),facecolor='none',edgecolors='C2')
# ax[1].scatter(lM200_T  ,(M200_z3 - M200_T)/np.mean(eM200_z3,axis=0),facecolor='none',edgecolors='C3')




# ax[1].plot([12.3,15],[0.,0.],'C7--')
# ax[1].axis([12.3,14.9,-2.5,2.5])

ax[0].set_xlabel(r'$\log (\langle M_{WL} \rangle)$')
# ax[1].set_ylabel(r'$(M^{ind}_{WL} - M_{WL})/ \epsilon(M^{ind}_{WL})$')

plt.savefig(path+'Rgroups/plots_newanalysis/M200_ztests.pdf',bbox_inches='tight')

# -----------------------
# -----------------------
# TESTCATS

f, ax = plt.subplots(figsize=(6.2,5))


ax.scatter(pcc_T,pcc_o4,facecolor='none',edgecolors='C0'  ,label='ODDS $> 0.4$')
ax.scatter(pcc_T,pcc_o6,facecolor='none',edgecolors='C1'  ,label='ODDS $> 0.6$')
ax.scatter(pcc_T,pcc_z2,facecolor='none',edgecolors='C2'  ,label='$z_{cut} = 0.2$')
ax.scatter(pcc_T,pcc_z3,facecolor='none',edgecolors='C3'  ,label='$z_{cut} = 0.3$')


ax.errorbar(pcc_T,pcc_o4,xerr=epcc_T,yerr=epcc_o4,facecolor='none',ecolor='C0'  ,fmt = 'none')
ax.errorbar(pcc_T,pcc_o6,xerr=epcc_T,yerr=epcc_o6,facecolor='none',ecolor='C1'  ,fmt = 'none')
ax.errorbar(pcc_T,pcc_z2,xerr=epcc_T,yerr=epcc_z2,facecolor='none',ecolor='C2'  ,fmt = 'none')
ax.errorbar(pcc_T,pcc_z3,xerr=epcc_T,yerr=epcc_z3,facecolor='none',ecolor='C3'  ,fmt = 'none')


ax.plot([0.2,1.0],[0.2,1.0],'C7--')
ax.legend(frameon = False,fontsize = 13)

ax.set_ylabel('$p^{ind}_{cc}$')
ax.set_xlabel('$p_{cc}$')

plt.savefig(path+'Rgroups/plots_newanalysis/pcc_ztests.pdf',bbox_inches='tight')

# -----------------------

# -----------------------
f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(lMH_T,lM200_T,facecolor='none',edgecolors='k'  ,label='Total sample')
ax[0].scatter(lMHc_T,lM200c_T,facecolor='none',edgecolors='C9'  ,label='C-sample')

ax[0].errorbar(lMH_T,lM200_T,yerr=elM200_T,fmt = 'none',ecolor='k')
ax[0].errorbar(lMHc_T,lM200c_T,yerr=elM200c_T,fmt = 'none',ecolor='C9')


ax[0].plot([12.3,15],[12.3,15],'C7--')
ax[0].legend(frameon = False,fontsize = 13)

ax[0].set_ylabel('$\log (M_{WL})$')
ax[0].axis([12.5,14.9,12.3,14.9])

ax[1].scatter(lMH_T  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
ax[1].scatter(lMHc_T  ,M200c_T/MHc_T,facecolor='none',edgecolors='C9',label='Total sample')

ax[1].errorbar(lMH_T  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1].errorbar(lMHc_T  ,M200c_T/MHc_T,yerr=eM200c_T/MHc_T,fmt = 'none',ecolor='C9')


ax[1].plot([12.3,15],[1.,1.],'C7--')
ax[1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')

plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_M200_Mbin_cM.pdf',bbox_inches='tight')

# -----------------------
plt.figure()

plt.scatter(lMH_T  ,ratio_cT,facecolor='none',edgecolors='k')
plt.scatter(lMH[mN1]  ,ratio_c[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMH[mN23] ,ratio_c[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMH[mN4M],ratio_c[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')
plt.scatter(lMH_H+0.05,ratio_cH,facecolor='none',edgecolors='C3',label='$  z \geq 0.13$')
plt.scatter(lMH_L+0.1,ratio_cL,facecolor='none',edgecolors='C0',label='$  z < 0.13$')

plt.errorbar(lMH_T  ,ratio_cT,yerr=eratio_cT,fmt = 'none',ecolor='k')
plt.errorbar(lMH[mN1]  ,ratio_c[mN1],yerr=eratio_c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMH[mN23] ,ratio_c[mN23],yerr=eratio_c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMH[mN4M],ratio_c[mN4M],yerr=eratio_c[:,mN4M],fmt = 'none',ecolor='C8')
plt.errorbar(lMH_H+0.05,ratio_cH,yerr=eratio_cH,fmt = 'none',ecolor='C3')
plt.errorbar(lMH_L+0.1,ratio_cL,yerr=eratio_cL,fmt = 'none',ecolor='C0')

plt.legend(frameon = False,fontsize = 14)
# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_c_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_c_t[j]+eratio_c_t[1][j],ratio_c_t[j]-eratio_c_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.5,14.8],[1.,1.],'C7--')
plt.axis([12.5,14.8,0,3.],'C7--')
plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$M^c_{WL}/M_{WL}$')
plt.savefig(path+'Rgroups/plots_newanalysis/ratioc.pdf',bbox_inches='tight')


# -----------------------hal
f, ax = plt.subplots(2, 2, figsize=(10,8),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0,0].scatter(lMH_T,lM200_T,facecolor='none',edgecolors='k',label='$N_{GAL} \geq 1$')
ax[0,0].scatter(lMH[mN1],lM200[mN1],facecolor='none',edgecolors='C6'  ,label='$N_{GAL} = 1$')
ax[0,0].scatter(lMH[mN23],lM200[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
ax[0,0].scatter(lMH[mN4M],lM200[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')
    
ax[0,0].errorbar(lMH_T,lM200_T,yerr=elM200_T,fmt = 'none',ecolor='k')
ax[0,0].errorbar(lMH[mN1],lM200[mN1],yerr=elM200[:,mN1],fmt = 'none',ecolor='C6'   )
ax[0,0].errorbar(lMH[mN23],lM200[mN23],yerr=elM200[:,mN23],fmt = 'none',ecolor='C7')
ax[0,0].errorbar(lMH[mN4M],lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8')
    
ax[0,0].plot([12.3,16],[12.3,16],'C7--')
ax[0,0].legend(loc =2,fontsize = 13)
    
ax[0,0].set_ylabel('$\log (M_{WL})$')
ax[0,0].axis([12.5,14.9,12.3,14.9])

ax[1,0].scatter(100.,100.,facecolor='none',edgecolors='k',label='Total sample')
ax[1,0].legend(fontsize = 13)

ax[1,0].scatter(lMH_T  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
ax[1,0].scatter(lMH[mN1]  ,(M200/MH)[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
ax[1,0].scatter(lMH[mN23] ,(M200/MH)[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
ax[1,0].scatter(lMH[mN4M],(M200/MH)[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')
    
ax[1,0].errorbar(lMH_T  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1,0].errorbar(lMH[mN1]  ,(M200/MH)[mN1],yerr=eM200[:,mN1]/MH[mN1],fmt = 'none',ecolor='C6')
ax[1,0].errorbar(lMH[mN23] ,(M200/MH)[mN23],yerr=eM200[:,mN23]/MH[mN23],fmt = 'none',ecolor='C7')
ax[1,0].errorbar(lMH[mN4M],(M200/MH)[mN4M],yerr=eM200[:,mN4M]/MH[mN4M],fmt = 'none',ecolor='C8')

ax[1,0].axis([12.5,14.9,0.2,1.6])
    
ax[1,0].plot([12.3,15.5],[1.,1.],'C7--')
ax[1,0].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1,0].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')


ax[0,1].plot(lMHc_T,lM200c_T        ,'ko'  ,label='Total sample')
ax[0,1].plot(lMHc[mN1],lM200c[mN1]  ,'C6o'  ,label='$N_{GAL} = 1$')
ax[0,1].plot(lMHc[mN23],lM200c[mN23],'C7o',label='$ 2 \leq N_{GAL} \leq 3$')
ax[0,1].plot(lMHc[mN4M],lM200c[mN4M],'C8o',label='$  N_{GAL} \geq 4$')
    
ax[0,1].errorbar(lMHc_T,lM200c_T,yerr=elM200c_T,fmt = 'none',ecolor='k')
ax[0,1].errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6'   )
ax[0,1].errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
ax[0,1].errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')
    
ax[0,1].plot([12.3,16],[12.3,16],'C7--')
ax[0,1].axis([12.5,14.9,12.3,14.9])
ax[0,1].set_yticks([])
   
ax[1,1].plot(100.,100.,'ko',label='$C-$sample')
ax[1,1].legend(fontsize = 13)
   
   
ax[1,1].plot(lMHc_T  ,M200c_T/MHc_T       ,'ko',label='Total sample')
ax[1,1].plot(lMHc[mN1]  ,(M200c/MHc)[mN1] ,'C6o',label='$N_{GAL} = 1$')
ax[1,1].plot(lMHc[mN23] ,(M200c/MHc)[mN23],'C7o',label='$ 2 \leq N_{GAL} \leq 3$')
ax[1,1].plot(lMHc[mN4M],(M200c/MHc)[mN4M] ,'C8o',label='$  N_{GAL} \geq 4$')
   
ax[1,1].errorbar(lMHc_T  ,M200c_T/MHc_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1,1].errorbar(lMHc[mN1]  ,(M200c/MHc)[mN1],yerr=eM200[:,mN1]/MH[mN1],fmt = 'none',ecolor='C6')
ax[1,1].errorbar(lMHc[mN23] ,(M200c/MHc)[mN23],yerr=eM200[:,mN23]/MH[mN23],fmt = 'none',ecolor='C7')
ax[1,1].errorbar(lMHc[mN4M],(M200c/MHc)[mN4M],yerr=eM200[:,mN4M]/MH[mN4M],fmt = 'none',ecolor='C8')

ax[1,1].axis([12.5,14.9,0.2,1.6])   
ax[1,1].plot([12.3,15.5],[1.,1.],'C7--')
ax[1,1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
# ax[1,1].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')
ax[1,1].set_yticks([])

ax[0,0].axvline(13.3,alpha=0.3,c='k')
ax[0,0].axvline(14.2,alpha=0.3,c='k')
ax[0,1].axvline(13.3,alpha=0.3,c='k')
ax[0,1].axvline(14.2,alpha=0.3,c='k')
ax[1,0].axvline(13.3,alpha=0.3,c='k')
ax[1,0].axvline(14.2,alpha=0.3,c='k')
ax[1,1].axvline(13.3,alpha=0.3,c='k')
ax[1,1].axvline(14.2,alpha=0.3,c='k')

plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_M200_Nbin.pdf',bbox_inches='tight')

# -----------------------
f, ax = plt.subplots(2, 2, figsize=(10,8),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0,0].scatter(lMH_T,lM200_T,facecolor='none',edgecolors='k',label='$N_{GAL} \geq 1$')
ax[0,0].scatter(lMH_L,lM200_L,facecolor='none',edgecolors='C0',label='$ z < 0.13$')
ax[0,0].scatter(lMH_H,lM200_H,facecolor='none',edgecolors='C3'  ,label='$z \geq 0.13$')

ax[0,0].errorbar(lMH_T,lM200_T,yerr=elM200_T,fmt = 'none',ecolor='k')
ax[0,0].errorbar(lMH_L,lM200_L,yerr=elM200_L,fmt = 'none',ecolor='C0')
ax[0,0].errorbar(lMH_H,lM200_H,yerr=elM200_H,fmt = 'none',ecolor='C3')

# extras

# ax[0,0].plot(lMH_e[0],lM200_e[0],'rv')
# ax[0,0].plot(lMH_e[1],lM200_e[1],'bv')
# ax[0,0].plot(lMH_e[2],lM200_e[2],'r^')
# ax[0,0].plot(lMH_e[3],lM200_e[3],'b^')

# ax[0,0].errorbar(lMH_e[0],lM200_e[0],yerr=(elM200_e[0,0]),fmt = 'none',ecolor='r')
# ax[0,0].errorbar(lMH_e[1],lM200_e[1],yerr=(elM200_e[0,1]),fmt = 'none',ecolor='b')
# ax[0,0].errorbar(lMH_e[2],lM200_e[2],yerr=(elM200_e[0,2]),fmt = 'none',ecolor='r')
# ax[0,0].errorbar(lMH_e[3],lM200_e[3],yerr=(elM200_e[0,3]),fmt = 'none',ecolor='b')

#-----------------

ax[0,0].plot([12.3,16],[12.3,16],'C7--')
ax[0,0].legend(fontsize = 13)
    
ax[0,0].set_ylabel('$\log (M_{WL})$')
ax[0,0].axis([12.5,14.9,12.3,14.9])

ax[1,0].scatter(100.,100.,facecolor='none',edgecolors='k',label='Total sample')
ax[1,0].legend(fontsize = 13)

ax[1,0].scatter(lMH_T  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
ax[1,0].scatter(lMH_L  ,M200_L/MH_L,facecolor='none',edgecolors='C0',label='$z < 0.13$')
ax[1,0].scatter(lMH_H  ,M200_H/MH_H,facecolor='none',edgecolors='C3',label='$z \geq 0.13$')
    
ax[1,0].errorbar(lMH_T  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1,0].errorbar(lMH_L  ,M200_L/MH_L,yerr=eM200_L/MH_L,fmt = 'none',ecolor='C0')
ax[1,0].errorbar(lMH_H  ,M200_H/MH_H,yerr=eM200_H/MH_H,fmt = 'none',ecolor='C3')

# extras

# ax[1,0].plot(lMH_e[0],(M200_e/MH_e)[0],'rv')
# ax[1,0].plot(lMH_e[1],(M200_e/MH_e)[1],'bv')
# ax[1,0].plot(lMH_e[2],(M200_e/MH_e)[2],'r^')
# ax[1,0].plot(lMH_e[3],(M200_e/MH_e)[3],'b^')
   
# ax[1,0].errorbar(lMH_e[0],(M200_e/MH_e)[0],yerr=(eM200_e/MH_e).T[0,0],fmt = 'none',ecolor='r')
# ax[1,0].errorbar(lMH_e[1],(M200_e/MH_e)[1],yerr=(eM200_e/MH_e).T[1,0],fmt = 'none',ecolor='b')
# ax[1,0].errorbar(lMH_e[2],(M200_e/MH_e)[2],yerr=(eM200_e/MH_e).T[2,0],fmt = 'none',ecolor='r')
# ax[1,0].errorbar(lMH_e[3],(M200_e/MH_e)[3],yerr=(eM200_e/MH_e).T[3,0],fmt = 'none',ecolor='b')

#-----------------


ax[1,0].axis([12.5,14.9,0.2,1.6])
    
ax[1,0].plot([12.3,15.5],[1.,1.],'C7--')
ax[1,0].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1,0].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')


ax[0,1].plot(lMHc_T,lM200c_T,'ko')
ax[0,1].plot(lMH_Lc,lM200_LcM,'C0o')
ax[0,1].plot(lMH_Hc,lM200_HcM,'C3o')
     
ax[0,1].errorbar(lMHc_T,lM200c_T,yerr=elM200_T,fmt = 'none',ecolor='k')
ax[0,1].errorbar(lMH_Lc,lM200_LcM,yerr=elM200_L,fmt = 'none',ecolor='C0')
ax[0,1].errorbar(lMH_Hc,lM200_HcM,yerr=elM200_H,fmt = 'none',ecolor='C3')
    
ax[0,1].plot([12.3,16],[12.3,16],'C7--')
# ax[0,1].legend(frameon = False,fontsize = 13)
    
# ax[0,1].set_ylabel('$\log (M_{WL})$')
ax[0,1].axis([12.5,14.9,12.3,14.9])
ax[0,1].set_yticks([])
   
ax[1,1].plot(100.,100.,'ko',label='$C-$sample')
ax[1,1].legend(fontsize = 13)
   
   
ax[1,1].plot(lMHc_T  ,M200c_T/MHc_T  ,'ko',label='Total sample')
ax[1,1].plot(lMH_Lc  ,M200_LcM/MH_LcM,'C0o',label='$z < 0.13$')
ax[1,1].plot(lMH_Hc  ,M200_HcM/MH_HcM,'C3o',label='$z \geq 0.13$')
     
ax[1,1].errorbar(lMHc_T  ,M200c_T/MHc_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
ax[1,1].errorbar(lMH_Lc  ,M200_LcM/MH_LcM,yerr=eM200_L/MH_L,fmt = 'none',ecolor='C0')
ax[1,1].errorbar(lMH_Hc  ,M200_HcM/MH_HcM,yerr=eM200_H/MH_H,fmt = 'none',ecolor='C3')

ax[1,1].axis([12.5,14.9,0.2,1.6])   
ax[1,1].plot([12.3,15.5],[1.,1.],'C7--')
ax[1,1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
# ax[1,1].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')
ax[1,1].set_yticks([])

ax[0,0].axvline(13.3,alpha=0.3,c='k')
ax[0,0].axvline(14.2,alpha=0.3,c='k')
ax[0,1].axvline(13.3,alpha=0.3,c='k')
ax[0,1].axvline(14.2,alpha=0.3,c='k')
ax[1,0].axvline(13.3,alpha=0.3,c='k')
ax[1,0].axvline(14.2,alpha=0.3,c='k')
ax[1,1].axvline(13.3,alpha=0.3,c='k')
ax[1,1].axvline(14.2,alpha=0.3,c='k')


plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_M200_zHL.pdf',bbox_inches='tight')

#----------------------------

def y(x,alpha):
    return np.log10(1.e14*((x/500.)**(alpha)))
    
def y_more(x):
    
    Mh = np.array([])
    
    for disp in x:
        a = 2.11 - (disp/200.)
        b = 0.19
        c = 0.20
        # a = 1.99 - (disp/200.)
        # b = 0.49
        # c = -0.20
        
        f0 = 0.45
        alpha_f = 0.61
        
        lMlumM = (-b + np.sqrt(b**2 - 4.*a*c))/(2.*c)
        lMlumm = (-b - np.sqrt(b**2 - 4.*a*c))/(2.*c)
        lMlum  = np.log10(10**(max(lMlumM,lMlumm))* 1e10)
        
        Mh = np.append(Mh,(((lMlum - 10.63)/0.29) + 11.94))
        # Mh = np.append(Mh,(((lMlum -  9.42)/0.98) + 11.29))
        
    return Mh

x = np.arange(250,800,20)
fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax2 = ax1.twiny()
# ax1.plot(out[-1][mN4M],lM200[mN4M],'C8s',ms=8)
ax1.plot(sdisp,lM200[mN4M],'C8s',ms=8)
ax1.errorbar(sdisp,lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8')
# ax1.errorbar(out[-1][mN4M],lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8')
ax1.set_xlim(250, 700)
# ax2.set_xlim(250, 560)
# ax2.set_xticks(sdisp)
# ax2.set_xticklabels(['13.56', '13.96', '14.22'])
# ax2.set_xlabel('$\log (M_{vir}/(h^{-1}M_\odot))$')
ax1.set_ylabel('$\log (M_{WL}/(h^{-1}M_\odot))$')
ax1.set_xlabel('$\sigma_V$ [km s$^{-1}]$')
ax1.plot(x,y(x,2.09),'C0--')
# ax1.plot(x,y_more(x),'C2--')
ax1.fill_between(x,y(x,1.89)+0.2,y(x,1.89)-0.2,alpha=0.3,color='C0')
ax1.plot(x,y(x,3.),'C1')
plt.savefig(path+'Rgroups/plots_newanalysis/Mdyn_M200.pdf',bbox_inches='tight')
# -----------------------

plt.figure()

plt.plot(lMdyn,lM200[mN4M],'C8s',ms=8)
# plt.errorbar(lMdyn,lM200[mN4M],xerr=slMdyn,yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8',label='$  N_{GAL} \geq 4$')
plt.errorbar(lMdyn,lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8',label='$  N_{GAL} \geq 4$')

# plt.plot(lMdyn_c2,lM200c[mN4M],'C8^')
# plt.errorbar(lMdyn_c2,lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')

plt.xlabel('$\log (M_{dyn})$')
plt.ylabel('$\log (M_{WL})$')
plt.axis([13.5,14.35,12.8,14.4])
plt.savefig(path+'Rgroups/plots_newanalysis/Mdyn_M200_NMbin.pdf',bbox_inches='tight')



# -----------------------
plt.figure()

plt.scatter(lMH_T,pcc_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lMH[mN1],pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMH[mN23],pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMH[mN4M],pcc[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lMH_T,pcc_T,yerr=epcc_T,fmt = 'none',ecolor='k')
plt.errorbar(lMH[mN1],pcc[mN1],yerr=epcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMH[mN23],pcc[mN23],yerr=epcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMH[mN4M],pcc[mN4M],yerr=epcc[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot(lMH_simu,pcc_simu,'C3')
plt.axis([12.5,14.8,0.1,1.3])
plt.legend(frameon = False,ncol = 2)

plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$p_{cc}$')
plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_pcc.pdf',bbox_inches='tight')

# -----------------------

# -----------------------
plt.figure()

plt.scatter(lMH_T,pcc_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lMH,pcc,facecolor='none',edgecolors='k')
plt.scatter(lMH_H,pcc_H,facecolor='none',edgecolors='k')
plt.scatter(lMH_L,pcc_L,facecolor='none',edgecolors='k')

plt.errorbar(lMH_T,pcc_T,yerr=epcc_T,fmt = 'none',ecolor='k')
plt.errorbar(lMH,pcc,yerr=epcc,fmt = 'none',ecolor='k')
plt.errorbar(lMH_H,pcc_H,yerr=epcc_H,fmt = 'none',ecolor='k')
plt.errorbar(lMH_L,pcc_L,yerr=epcc_L,fmt = 'none',ecolor='k')

plt.scatter(lMHc_T+0.1,pcc_cT,facecolor='none',edgecolors='C9',label='C-sample')
plt.scatter(lMHc+0.1,pcc_c,facecolor='none',edgecolors='C9')
plt.scatter(lMH_H+0.1,pcc_Hc,facecolor='none',edgecolors='C9')
plt.scatter(lMH_L+0.1,pcc_Lc,facecolor='none',edgecolors='C9')

plt.errorbar(lMHc_T+0.1,pcc_cT,yerr=epcc_cT,fmt = 'none',ecolor='C9')
plt.errorbar(lMHc+0.1,pcc_c,yerr=epcc_c,fmt = 'none',ecolor='C9')
plt.errorbar(lMH_H+0.1,pcc_Hc,yerr=epcc_Hc,fmt = 'none',ecolor='C9')
plt.errorbar(lMH_L+0.1,pcc_Lc,yerr=epcc_Lc,fmt = 'none',ecolor='C9')


plt.plot(lMH_simu,pcc_simu,'C3')
plt.axis([12.6,15,0.,1.5])
plt.legend(frameon = False,ncol = 2)

plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$p_{cc}$')
plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_pcc_cM.pdf',bbox_inches='tight')

# -----------------------


plt.figure()
plt.hist(np.append(out[14],out_T[14]),7,color='k',histtype='step')
plt.hist(np.append(outc[14],outc_T[14]),7,color='C9',histtype='step',label='$C-$sample')
plt.legend(fontsize=14,loc=2,frameon=False)
plt.ylabel('N')
plt.xlabel('$p_{cc}$')
plt.savefig(path+'Rgroups/plots_newanalysis/pcc_dist.pdf',bbox_inches='tight')


# -----------------------
plt.figure()


plt.scatter(out_T[9]  ,M200_T/MH_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(out[9][mN1]  ,(M200/MH)[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(out[9][mN23] ,(M200/MH)[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(out[9][mN4M],(M200/MH)[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')
plt.scatter(out_L[9]  ,M200_L/MH_L,facecolor='none',edgecolors='C0',label='$z < 0.13$')
plt.scatter(out_H[9]  ,M200_H/MH_H,facecolor='none',edgecolors='C3',label='$z \geq 0.13$')

plt.errorbar(out_T[9]  ,M200_T/MH_T,yerr=eM200_T/MH_T,fmt = 'none',ecolor='k')
plt.errorbar(out[9][mN1]  ,(M200/MH)[mN1],yerr=eM200[:,mN1]/MH[mN1],fmt = 'none',ecolor='C6')
plt.errorbar(out[9][mN23] ,(M200/MH)[mN23],yerr=eM200[:,mN23]/MH[mN23],fmt = 'none',ecolor='C7')
plt.errorbar(out[9][mN4M],(M200/MH)[mN4M],yerr=eM200[:,mN4M]/MH[mN4M],fmt = 'none',ecolor='C8')
plt.errorbar(out_L[9]  ,M200_L/MH_L,yerr=eM200_L/MH_L,fmt = 'none',ecolor='C0')
plt.errorbar(out_H[9]  ,M200_H/MH_H,yerr=eM200_H/MH_H,fmt = 'none',ecolor='C3')

plt.legend(loc= 2,frameon = False,fontsize = 13,ncol=2)

plt.plot([0.06,0.175],[1.,1.],'C7--')
plt.xlabel(r'$\langle z \rangle$')
plt.ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')
plt.savefig(path+'Rgroups/plots_newanalysis/ratio_MH_z.pdf',bbox_inches='tight')


# -----------------------
plt.figure()


plt.plot(outc_T[9]  ,M200c_T/MHc_T       ,'ko',label='Total sample')
plt.plot(outc[9][mN1]  ,(M200c/MHc)[mN1] ,'C6o',label='$N_{GAL} = 1$')
plt.plot(outc[9][mN23] ,(M200c/MHc)[mN23],'C7o',label='$ 2 \leq N_{GAL} \leq 3$')
plt.plot(outc[9][mN4M],(M200c/MHc)[mN4M] ,'C8o',label='$  N_{GAL} \geq 4$')
plt.plot(out_LcM[9]  ,M200_LcM/MH_LcM,'C0o',label='$z < 0.13$')
plt.plot(out_HcM[9]  ,M200_HcM/MH_HcM,'C3o',label='$z \geq 0.13$')

plt.errorbar(outc_T[9]  ,M200c_T/MHc_T,yerr=eM200c_T/MHc_T,fmt = 'none',ecolor='k')
plt.errorbar(outc[9][mN1]  ,(M200c/MHc)[mN1],yerr=eM200c[:,mN1]/MHc[mN1],fmt = 'none',ecolor='C6')
plt.errorbar(outc[9][mN23] ,(M200c/MHc)[mN23],yerr=eM200c[:,mN23]/MHc[mN23],fmt = 'none',ecolor='C7')
plt.errorbar(outc[9][mN4M],(M200c/MHc)[mN4M],yerr=eM200c[:,mN4M]/MHc[mN4M],fmt = 'none',ecolor='C8')
plt.errorbar(out_LcM[9]  ,M200_LcM/MH_LcM,yerr=eM200_LcM/MH_LcM,fmt = 'none',ecolor='C0')
plt.errorbar(out_HcM[9]  ,M200_HcM/MH_HcM,yerr=eM200_HcM/MH_HcM,fmt = 'none',ecolor='C3')


plt.legend(loc= 2,frameon = False,fontsize = 13,ncol=2)

plt.plot([0.06,0.175],[1.,1.],'C7--')
plt.xlabel(r'$\langle z \rangle$')
plt.ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')
plt.savefig(path+'Rgroups/plots_newanalysis/ratio_MH_z_cM.pdf',bbox_inches='tight')


'''
# -----------------------

# -----------------------
f, ax = plt.subplots(2, 1, figsize=(6.2,10),sharex=True)
f.subplots_adjust(hspace=0,wspace=0)


ax[0].plot(lMH[mN1],lM200[mN1]  ,'ko'  ,label='$N_{GAL} = 1$')
ax[0].plot(lMH[mN23],lM200[mN23],'ks',label='$ 2 \leq N_{GAL} \leq 3$')
ax[0].plot(lMH[mN4M],lM200[mN4M],'k^',label='$  N_{GAL} \geq 4$')

ax[0].errorbar(lMH[mN1],lM200[mN1],yerr=elM200[:,mN1],fmt = 'none',ecolor='k'   )
ax[0].errorbar(lMH[mN23],lM200[mN23],yerr=elM200[:,mN23],fmt = 'none',ecolor='k')
ax[0].errorbar(lMH[mN4M],lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='k')

ax[0].legend(frameon = False,fontsize = 13)

ax[0].plot(lMHc[mN1] ,lM200c[mN1] ,'C9o'  ,label='$N_{GAL} = 1$')
ax[0].plot(lMHc[mN23],lM200c[mN23],'C9s',label='$ 2 \leq N_{GAL} \leq 3$')
ax[0].plot(lMHc[mN4M],lM200c[mN4M],'C9^',label='$  N_{GAL} \geq 4$')

ax[0].errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C9'   )
ax[0].errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C9')
ax[0].errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C9')


ax[0].plot([12.3,16],[12.3,16],'C7--')


ax[0].set_ylabel('$\log (M_{WL})$')
ax[0].axis([12.3,15.0,12.3,15.])


ax[1].plot(lMH[mN1]  ,(M200/MH)[mN1] ,'ko',label='$N_{GAL} = 1$')
ax[1].plot(lMH[mN23] ,(M200/MH)[mN23],'ks',label='$ 2 \leq N_{GAL} \leq 3$')
ax[1].plot(lMH[mN4M],(M200/MH)[mN4M] ,'k^',label='$  N_{GAL} \geq 4$')

ax[1].errorbar(lMH[mN1]  ,(M200/MH)[mN1],yerr=eM200[:,mN1]/MH[mN1],fmt = 'none',ecolor='k')
ax[1].errorbar(lMH[mN23] ,(M200/MH)[mN23],yerr=eM200[:,mN23]/MH[mN23],fmt = 'none',ecolor='k')
ax[1].errorbar(lMH[mN4M],(M200/MH)[mN4M],yerr=eM200[:,mN4M]/MH[mN4M],fmt = 'none',ecolor='k')

ax[1].plot(lMHc[mN1]  ,(M200c/MHc)[mN1] ,'C9o',label='$N_{GAL} = 1$')
ax[1].plot(lMHc[mN23] ,(M200c/MHc)[mN23],'C9s',label='$ 2 \leq N_{GAL} \leq 3$')
ax[1].plot(lMHc[mN4M],(M200c/MHc)[mN4M] ,'C9^',label='$  N_{GAL} \geq 4$')

ax[1].errorbar(lMHc[mN1]  ,(M200c/MHc)[mN1],yerr=eM200c[:,mN1]/MHc[mN1],fmt = 'none',ecolor='C9')
ax[1].errorbar(lMHc[mN23] ,(M200c/MHc)[mN23],yerr=eM200c[:,mN23]/MHc[mN23],fmt = 'none',ecolor='C9')
ax[1].errorbar(lMHc[mN4M],(M200c/MHc)[mN4M],yerr=eM200c[:,mN4M]/MHc[mN4M],fmt = 'none',ecolor='C9')


ax[1].plot([12.3,15.5],[1.,1.],'C7--')
ax[1].set_xlabel(r'$\log (\langle M_{AM} \rangle)$')
ax[1].set_ylabel(r'$M_{WL}/\langle M_{AM} \rangle$')

plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_M200_NMbin_cM.pdf',bbox_inches='tight')


plt.figure()
# plt.plot(100,100,'k^', label = '$C-$sample')
# plt.scatter(100,100,facecolor='none',edgecolors='k', label = 'Total sample')

plt.scatter(10.*out_T[6]    ,10.*out_T[11],facecolor='none',edgecolors='k'  ,label='Total sample')
plt.scatter(10.*out[6][mN1] ,10.*out[11][mN1],facecolor='none',edgecolors='C6'  ,label='$N_{GAL} = 1$')
plt.scatter(10.*out[6][mN23],10.*out[11][mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(10.*out[6][mN4M],10.*out[11][mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(10.*out_T[6]    ,10.*out_T[11],yerr=10.*eM200_T,fmt = 'none',ecolor='k')
plt.errorbar(10.*out[6][mN1] ,10.*out[11][mN1],yerr=10.*eM200[:,mN1],fmt = 'none',ecolor='C6'  )
plt.errorbar(10.*out[6][mN23],10.*out[11][mN23],yerr=10.*eM200[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(10.*out[6][mN4M],10.*out[11][mN4M],yerr=10.*eM200[:,mN4M],fmt = 'none',ecolor='C8')

# plt.plot(lMHc[mN1] ,lM200c[mN1] ,'C6^')
# plt.plot(lMHc[mN23],lM200c[mN23],'C7^')
# plt.plot(lMHc[mN4M],lM200c[mN4M],'C8^')

# plt.errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6')
# plt.errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
# plt.errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([0,40.],[0,40],'C7--')
plt.legend(frameon = False,fontsize = 14)
plt.xlabel(r'$\langle M_{AM} \rangle [10^{13}M_\odot h^{-1}]$')
plt.ylabel('$M_{WL} [10^{13}M_\odot h^{-1}]$')
plt.axis([0,35,0,40])
plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_M200_NMbin_lin.pdf',bbox_inches='tight')

# -----------------------
# -----------------------
plt.figure()

plt.scatter(lM200_T,pcc_T,facecolor='none',edgecolors='k',label='Total sample')
plt.scatter(lM200[mN1],pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23],pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],pcc[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200_T,pcc_T,xerr=elM200_T,yerr=epcc_T,fmt = 'none',ecolor='k')
plt.errorbar(lM200[mN1],pcc[mN1],xerr=elM200[:,mN1],yerr=epcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23],pcc[mN23],xerr=elM200[:,mN23],yerr=epcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],pcc[mN4M],xerr=elM200[:,mN4M],yerr=epcc[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot(lMH_simu,pcc_simu,'C3')
plt.axis([12.5,14.5,0.1,1.2])
plt.legend(frameon = False,ncol = 2)
# plt.plot(lM200c[mN1] ,pcc_c[mN1]  ,'C6^',label='$N_{GAL} = 1$')
# plt.plot(lM200c[mN23],pcc_c[mN23],'C7^',label='$ 2 \leq N_{GAL} \leq 3$')
# plt.plot(lM200c[mN4M],pcc_c[mN4M],'C8^',label='$ 4 \leq N_{GAL} \leq 10$')

# plt.errorbar(lM200c[mN1] ,pcc_c[mN1] ,xerr=elM200c[:,mN1] ,yerr=epcc_c[:,mN1],fmt = 'none',ecolor='C6')
# plt.errorbar(lM200c[mN23],pcc_c[mN23],xerr=elM200c[:,mN23],yerr=epcc_c[:,mN23],fmt = 'none',ecolor='C7')
# plt.errorbar(lM200c[mN4M],pcc_c[mN4M],xerr=elM200c[:,mN4M],yerr=epcc_c[:,mN4M],fmt = 'none',ecolor='C8')


# plt.legend(frameon = False)
plt.xlabel('$\log (M_{WL})$')
plt.ylabel('$p_{cc}$')
plt.savefig(path+'Rgroups/plots_newanalysis/M200_pcc_NMbin.pdf',bbox_inches='tight')

# -----------------------





plt.figure()

plt.scatter(lM200[mN1]  ,ratio_pcc[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23] ,ratio_pcc[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],ratio_pcc[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200[mN1]  ,ratio_pcc[mN1],xerr=elM200[:,mN1],yerr=eratio_pcc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23] ,ratio_pcc[mN23],xerr=elM200[:,mN23],yerr=eratio_pcc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],ratio_pcc[mN4M],xerr=elM200[:,mN4M],yerr=eratio_pcc[:,mN4M],fmt = 'none',ecolor='C8')

# j = np.argsort(lM200c_t)
# plt.plot(lM200c_t[j],ratio_pcc_t[j],'C4')
# plt.fill_between(lM200c_t[j],ratio_pcc_t[j]+eratio_pcc_t[1][j],ratio_pcc_t[j]-eratio_pcc_t[0][j],color='C4',alpha=0.5)


# plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.3],[1.,1.],'C7--')
plt.xlabel('$\log (M_{WL})$')
plt.ylabel('$p^c_{cc}/p_{cc}$')
plt.savefig(path+'Rgroups/plots_newanalysis/ratiopcc_NMbin.pdf',bbox_inches='tight')

# -----------------------





# -----------------------
plt.figure()
plt.plot(lMHc_t,lM200c_t,'C4')
plt.scatter(lMHc[mN1],lM200c[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lMHc[mN23],lM200c[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lMHc[mN4M],lM200c[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.fill_between(lMHc_t,lM200c_t-elM200c_t[0],lM200_t+elM200c_t[1],color='C4',alpha=0.5)
# plt.errorbar(lMHc_t,lM200c_t,yerr=elM200c_t,fmt = 'none',ecolor='C0')
plt.errorbar(lMHc[mN1],lM200c[mN1],yerr=elM200c[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lMHc[mN23],lM200c[mN23],yerr=elM200c[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lMHc[mN4M],lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot([12.3,14.6],[12.3,14.6],'C7--')
plt.legend(frameon = False)
plt.xlabel('$\log (M_{AM})$')
plt.ylabel('$\log (M_{WL})$')
plt.axis([12.3,14.2,12.3,14.2])
plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_M200_NMbin_conM.pdf',bbox_inches='tight')
# -----------------------

# -----------------------
plt.figure()

plt.scatter(lM200[mN1]  ,ratio_LH[mN1],facecolor='none',edgecolors='C6',label='$N_{GAL} = 1$')
plt.scatter(lM200[mN23] ,ratio_LH[mN23],facecolor='none',edgecolors='C7',label='$ 2 \leq N_{GAL} \leq 3$')
plt.scatter(lM200[mN4M],ratio_LH[mN4M],facecolor='none',edgecolors='C8',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200[mN1]  ,ratio_LH[mN1],xerr=elM200[:,mN1],yerr=eratio_LH[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200[mN23] ,ratio_LH[mN23],xerr=elM200[:,mN23],yerr=eratio_LH[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200[mN4M],ratio_LH[mN4M],xerr=elM200[:,mN4M],yerr=eratio_LH[:,mN4M],fmt = 'none',ecolor='C8')

plt.plot(lM200c[mN1]  ,ratio_LHc[mN1],'C6^',label='$N_{GAL} = 1$')
plt.plot(lM200c[mN23] ,ratio_LHc[mN23],'C7^',label='$ 2 \leq N_{GAL} \leq 3$')
plt.plot(lM200c[mN4M],ratio_LHc[mN4M],'C8^',label='$  N_{GAL} \geq 4$')

plt.errorbar(lM200c[mN1]  ,ratio_LHc[mN1],xerr=elM200c[:,mN1],yerr=eratio_LHc[:,mN1],fmt = 'none',ecolor='C6')
plt.errorbar(lM200c[mN23] ,ratio_LHc[mN23],xerr=elM200c[:,mN23],yerr=eratio_LHc[:,mN23],fmt = 'none',ecolor='C7')
plt.errorbar(lM200c[mN4M],ratio_LHc[mN4M],xerr=elM200c[:,mN4M],yerr=eratio_LHc[:,mN4M],fmt = 'none',ecolor='C8')

j = np.argsort(lM200c_t)
plt.plot(lM200c_t[j],ratio_LHc_t[j],'C4')
plt.fill_between(lM200c_t[j],ratio_LHc_t[j]+eratio_LHc_t[1][j],ratio_LHc_t[j]-eratio_LHc_t[0][j],color='C4',alpha=0.5)

# plt.plot(lM200_t,ratio_LH_t,'k')
# plt.fill_between(lM200_t,ratio_LH_t-eratio_LH_t[0],ratio_LH_t+eratio_LH_t[1],color='k',alpha=0.5)


plt.yscale('log')

# plt.legend(frameon = False)
plt.plot([12.3,14.3],[1.,1.],'C7--')
plt.xlabel('$\log (M_{WL})$')
plt.ylabel('$M^H_{200}/M^L_{200}$')
plt.savefig(path+'Rgroups/plots_newanalysis/ratioLH_NMbin.pdf',bbox_inches='tight')

# -----------------------
# -----------------------
plt.figure()

plt.scatter(lMdyn,lM200[mN4M],facecolor='none',edgecolors='C8')
plt.errorbar(lMdyn,lM200[mN4M],yerr=elM200[:,mN4M],fmt = 'none',ecolor='C8',label='$  N_{GAL} \geq 4$')

plt.plot(lMdyn_c,lM200c[mN4M],'C8^')
plt.errorbar(lMdyn_c,lM200c[mN4M],yerr=elM200c[:,mN4M],fmt = 'none',ecolor='C8')


plt.plot([12.3,14.6],[12.3,14.6],'C7--')

plt.xlabel('$\log (M_{dyn})$')
plt.ylabel('$\log (M_{WL})$')
plt.axis([13.,14.5,12.7,14.5])
plt.savefig(path+'Rgroups/plots_newanalysis/Mdyn_M200_NMbin.pdf',bbox_inches='tight')




# -----------------------
plt.figure()

plt.scatter(lMH_H,pcc_H,facecolor='none',edgecolors='C3',label='High z')
plt.scatter(lMH_L,pcc_L,facecolor='none',edgecolors='C0',label='Low z')

plt.errorbar(lMH_H,pcc_H,yerr=epcc_H,fmt = 'none',ecolor='C3')
plt.errorbar(lMH_L,pcc_L,yerr=epcc_L,fmt = 'none',ecolor='C0')

plt.xlabel(r'$\log (\langle M_{AM} \rangle)$')
plt.ylabel('$p_{cc}$')
plt.savefig(path+'Rgroups/plots_newanalysis/Mhalo_pcc_zbin.pdf',bbox_inches='tight')
'''
