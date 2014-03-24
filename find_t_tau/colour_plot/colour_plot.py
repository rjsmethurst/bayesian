import numpy as N
import pylab as P
import scipy as S
import pyfits as F
from t_tau_func import *
from astropy.cosmology import FlatLambdaCDM

#cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
#av_z = 0.076401
#age = cosmo.age(av_z).value
age=13.7

font = {'family':'serif', 'size':16}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')

#dir ='/home/smethurst/Projects/Green-Valley-Project/bc03/models/Padova1994/chabrier/ASCII/'
#model = 'extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
#data = N.loadtxt(dir+model)
#
tq = N.linspace(0.0, 13.8, 200)
tau = N.linspace(0.001, 3, 200)

time = N.arange(0, 0.01, 0.003)
t = N.linspace(0, 13.7, 200)
t = N.append(time, t[1:])


#ur = N.zeros((len(tau),len(tq)))
#nuv = N.zeros_like(ur)
#for n in range(len(tq)):
#    for m in range(len(tau)):
#        nuv[m,n], ur[m,n] = predict_c_one([tq[n], tau[m]], age)
#
#N.save('ur.npy', ur)
#N.save('nuvu.npy', nuv)

ur = N.load('ur.npy')
nuv = N.load('nuvu.npy')

sfr = N.zeros((len(tau), len(tq)))
for n in range(len(tq)):
    for m in range(len(tau)):
        full_sfr = expsfh(tau[m], tq[n], t)
        sfr[m,n] = 100 - (N.interp(age, t, full_sfr)/full_sfr[0])*100

#P.figure()
#P.imshow(sfr, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
#P.text(0.5, 2.75, r'$t_{obs} = 13.7 ~Gyr$')
#P.xlabel(r'$t_{quench} (Gyr)$')
#P.ylabel(r'$\tau$ (Gyr)')
#cbar = P.colorbar()
#cbar.set_label(r'% drop in star formation rate')
#P.savefig('sfr_drop.pdf')

fig = P.figure(figsize=(15,4))
P.subplot(1,3,1)
ax1 = P.imshow(ur, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar = P.colorbar(ax1)
cbar.set_label(r'predicted $u-r$ colour')

P.subplot(1,3,2)
ax2 = P.imshow(nuv, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar= P.colorbar(ax2)
cbar.set_label(r'predicted $NUV-u$ colour')

P.subplot(1,3,3)
ax3 = P.imshow(sfr, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
P.text(0.5, 2.75, r'$t_{obs} = 13.7 ~Gyr$')
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar = P.colorbar(ax3)
cbar.set_label(r'% drop in SFR')

save_fig = '/home/smethurst/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_plot/sfr_drop_and_colours_'+str(len(tq))+'_t_tau_obs_'+str(age)+'_Gyr_with_203_timesteps.pdf'
fig.tight_layout()
fig.savefig(save_fig)

        