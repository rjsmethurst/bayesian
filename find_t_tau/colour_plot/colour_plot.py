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

fig = P.figure(figsize=(15,6))
P.subplot(1,3,1)
ax1 = P.imshow(ur, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], alpha = 0.7, cmap=P.cm.spectral)
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar = P.colorbar(ax1, orientation='horizontal')
cbar.set_label(r'predicted $u-r$ colour', labelpad=10)
cbar.set_ticks([0.75, 1.25, 1.75, 2.25])
cbar.set_ticklabels([0.75, 1.25, 1.75, 2.25])

P.subplot(1,3,2)
ax2 = P.imshow(nuv, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], alpha = 0.7, cmap=P.cm.spectral)
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar= P.colorbar(ax2, orientation='horizontal')
cbar.set_label(r'predicted $NUV-u$ colour', labelpad=10)
cbar.set_ticks([0.4, 1.2, 2.0, 2.8])
cbar.set_ticklabels([0.4, 1.2, 2.0, 2.8])


P.subplot(1,3,3)
ax3 = P.imshow(sfr, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], alpha = 0.7)
P.text(0.5, 2.75, r'$t_{obs} = 13.7 ~Gyr$')
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar = P.colorbar(ax3, orientation='horizontal')
cbar.set_label(r'% drop in SFR', labelpad=10)
cbar.set_ticks([0, 20, 40, 60, 80, 100])
cbar.set_ticklabels([0, 20, 40, 60, 80, 100])

save_fig = '/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_plot/sfr_drop_and_colours_'+str(len(tq))+'_t_tau_obs_'+str(age)+'_Gyr_with_203_timesteps.pdf'
fig.tight_layout()
fig.savefig(save_fig)

        