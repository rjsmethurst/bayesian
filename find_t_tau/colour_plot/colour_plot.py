import numpy as N
import pylab as P
import scipy as S
import pyfits as F
from t_tau_func import *
from astropy.cosmology import FlatLambdaCDM

#cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
#av_z = 0.076401
#age = cosmo.age(av_z).value
age=12.878505072906682

font = {'family':'serif', 'size':16}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')


#
tq = N.linspace(0.0, 13.8, 75)
tau = N.linspace(0.001, 3, 75)

tqs = N.outer(tq, N.ones(len(tq)))
taus = N.outer(N.ones(len(tau)), tau)

time = N.arange(0, 0.01, 0.003)
t = N.linspace(0, 13.7, 100)
t = N.append(time, t[1:])

#
#ur = N.zeros((len(tau),len(tq)))
#nuv = N.zeros_like(ur)
#for n in range(len(tq)):
#    for m in range(len(tau)):
#        nuv[m,n], ur[m,n] = predict_c_one([tq[n], tau[m]], age)
#
#N.save('ur_new_sfhs_avg_age.npy', ur)
#N.save('nuvu_new_sfhs_avg_age.npy', nuv)

ur = N.load('ur_new_sfhs_avg_age.npy')
nuv = N.load('nuvu_new_sfhs_avg_age.npy')

sfr = N.zeros((len(tau), len(tq)))
mass = N.zeros_like(sfr)
ssfr = N.zeros_like(sfr)
for n in range(len(tq)):
    for m in range(len(tau)):
        full_sfr, full_mass = expsfh(tau[m], tq[n], t)
        sfr[m,n] = N.interp(age, t, full_sfr)
        mass[m,n] = N.interp(age, t, full_mass)
        ssfr[m,n] = sfr[m,n]/mass[m,n]

N.save('sfr_new_sfhs_avg_age.npy', sfr)
N.save('mass_new_sfhs_avg_age.npy', mass)

P.figure()
P.scatter(N.log10(mass), N.log10(sfr), c=tqs, s=taus*10, linewidth=0.1)
P.ylim(ymin=-2.5, ymax = 2)
P.xlim((8.5,12.0))
P.colorbar()
P.show()

print sfr
print ssfr

#P.figure()
#P.imshow(sfr, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
#P.text(0.5, 2.75, r'$t_{obs} = 13.7 ~Gyr$')
#P.xlabel(r'$t_{quench} (Gyr)$')
#P.ylabel(r'$\tau$ (Gyr)')
#cbar = P.colorbar()
#cbar.set_label(r'% drop in star formation rate')
#P.savefig('sfr_drop.pdf')

fig = P.figure(figsize=(6.5,15))
ax1= P.subplot(3,1,1)
im = ax1.imshow(ur, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], alpha = 0.5, cmap=P.cm.spectral)
[l.set_rotation(45) for l in ax1.get_xticklabels()]
[j.set_rotation(45) for j in ax1.get_yticklabels()]
P.text(0.5, 2.75, r'$t_{obs} = 12.8  Gyr$')
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar = P.colorbar(im, orientation='vertical')
cbar.set_label(r'predicted $u-r$ colour', labelpad=10)
#cbar.set_ticks([0.75, 1.25, 1.75, 2.25])
#cbar.set_ticklabels([0.75, 1.25, 1.75, 2.25])

ax2 =P.subplot(3,1,2)
im = ax2.imshow(nuv, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], alpha = 0.5, cmap=P.cm.spectral)
[l.set_rotation(45) for l in ax2.get_xticklabels()]
[j.set_rotation(45) for j in ax2.get_yticklabels()]
P.text(0.5, 2.75, r'$t_{obs} = 12.8  Gyr$')
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar= P.colorbar(im, orientation='vertical')
cbar.set_label(r'predicted $NUV-u$ colour', labelpad=10)
#cbar.set_ticks([0.4, 1.2, 2.0, 2.8])
#cbar.set_ticklabels([0.4, 1.2, 2.0, 2.8])


ax3 = P.subplot(3,1,3)
im = ax3.imshow(sfr, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], alpha = 0.5, cmap=P.cm.spectral_r)
[l.set_rotation(45) for l in ax3.get_xticklabels()]
[j.set_rotation(45) for j in ax3.get_yticklabels()]
P.text(0.5, 2.75, r'$t_{obs} = 12.8  Gyr$')
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$ (Gyr)')
cbar = P.colorbar(im, orientation='vertical')
cbar.set_label(r'$ SFR [M_{\odot} yr^{-1}]$', labelpad=10)
#cbar.set_ticks([0.008, 0.024, 0.040, 0.056, 0.072])
#cbar.set_ticklabels([0.008, 0.024, 0.040, 0.056, 0.072])
save_fig = '/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_plot/new_sfhs_sfr_drop_and_colours_'+str(len(tq))+'_t_tau_obs_'+str(age)+'_Gyr_with_203_timesteps.pdf'
fig.tight_layout()
fig.savefig(save_fig)

        