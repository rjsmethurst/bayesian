import numpy as N
import pylab as P
import scipy as S
import pyfits as F
from t_tau_func import *
from astropy.cosmology import FlatLambdaCDM
import triangle

#cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
#av_z = 0.076401
#age = cosmo.age(av_z).value
age=N.array([13.8, 12.8, 10, 8, 6])

font = {'family':'serif', 'size':16}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')

filem = '/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match_GZ1_k_correct_MPA_JHU_SFR_mass.fits'
datm = F.open(filem)
gz2datam = datm[1].data
avg_mass = gz2datam.field('AVG_MASS')
avg_sfr = gz2datam.field('AVG_SFR')
u_r = gz2datam.field('MU_MR')
nuv_u = gz2datam.field('NUV_U')

time = N.arange(0, 0.01, 0.003)
t = N.linspace(0, 13.8, 100)
t = N.append(time, t[1:])

tq_r = N.random.rand(10000)
tq_ri = N.random.randint(0.0,14.0,10000)
tq = tq_r + tq_ri

P.figure()
P.hist(tq, bins=10)
P.show()

tau_r = N.random.rand(10000)
tau_ri = N.random.randint(0,3,10000)
tau = tau_r + tau_ri

print N.max(tau)

P.figure()
P.hist(tau, bins=10)
P.show()

#full_sfr = N.load('full_sfr.npy')
#full_mass = N.load('full_mass.npy')

#ur = N.zeros((len(tau), len(age)))
#nuv = N.zeros_like(ur)
#sfr = N.zeros(len(tau)*len(age)).reshape(len(tau), len(age))
#mass = N.zeros_like(sfr)
#ssfr = N.zeros_like(sfr)
#dsfr = N.zeros_like(sfr)
##full_sfr = N.zeros((len(tau), len(t)))
##full_mass = N.zeros((len(tau), len(t)))
#for m in range(len(tau)):
#    nuv[m,:], ur[m,:], full_sfr[m,:], full_mass[m,:] = predict_c_one([tq[m], tau[m]], age)
#    sfr[m, :] = N.interp(age, t, full_sfr[m,:])
#    mass[m,:] = N.interp(age, t, full_mass[m,:])
#    ssfr[m,:] = sfr[m,:]/mass[m,:]
#    dsfr[m,:] = full_sfr[m,0] - sfr[m,:]
#N.save('ur_new_sfhs_avg_age.npy', ur)
#N.save('nuvu_new_sfhs_avg_age.npy', nuv)


ur = N.load('ur_new_sfhs_avg_age.npy')
nuv = N.load('nuvu_new_sfhs_avg_age.npy')

#N.save('sfr_new_sfhs_avg_age.npy', sfr)
#N.save('mass_new_sfhs_avg_age.npy', mass)
#
#N.save('full_sfr.npy', full_sfr)
#N.save('full_mass.npy', full_mass)

#sfr = N.load('sfr_new_sfhs_avg_age.npy')
#mass = N.load('mass_new_sfhs_avg_age.npy')

P.figure(figsize=(30,5))
for n in range(len(age)):
    ax = P.subplot(1, len(age)+1, n+2, aspect='auto')
    triangle.hist2d(ur[:,n], nuv[:,n], ax=ax,bins=50, extent=([0.0, 3.0],[-1.0,5.0]))
    [l.set_rotation(45) for l in ax.get_xticklabels()]
    [l.set_rotation(45) for l in ax.get_yticklabels()]
    ax.set_ylabel(r'predicted $NUV-u$')
    ax.set_xlabel(r'predicted $u-r$')
    ax.text(0.25, 4.5, 't = '+str(age[n])+' Gyr')
    print n
ax5 = P.subplot(1,6,1, aspect='auto')
triangle.hist2d(u_r, nuv_u, ax=ax5, bins=50, extent=([0.0, 3.0],[-1.0,5.0]))
ax5.set_xlabel(r'$u-r$')
ax5.set_ylabel(r'$NUV-u$')
[l.set_rotation(45) for l in ax5.get_xticklabels()]
[l.set_rotation(45) for l in ax5.get_yticklabels()]
ax5.text(0.25, 4.5, 'GZ2 galaxies')
P.tight_layout()
P.savefig('pred_colour.pdf')
P.show()

x = N.arange(8, 12.25, 0.25)
m = (1.0--1)/(11.0-9)
y = m*x + (-1-(9*m))
uy = y+0.3
ly = y-0.3

#P.figure()
#for m in range(len(tau)):
#    P.plot(N.log10(full_mass[m,:]), N.log10(full_sfr[m,:]))
#P.plot(x,y, color='r')
#P.plot(x, uy, color='r', linestyle='dashed')
#P.plot(x, ly, color='r', linestyle='dashed')
#P.xlim((8, 12.25))
#P.ylim((-2.5, 2.0))
#P.show()

log_sfr = N.log10(sfr)

fig = P.figure(figsize=(24,5))
for n in range(len(age)):
    ax = P.subplot(1, len(age)+1, n+2, aspect='auto')
    #ax.scatter(N.log10(mass[:,n]), sfr[:,n], rasterized='True')
    ax.text(8.3, 1.5, 't = '+str(age[n])+' Gyr')
    triangle.hist2d(N.log10(mass[:,n]), log_sfr[:,n], ax=ax, bins=50, extent=([8.0, 12.0],[-2.5, N.max(log_sfr[:,n])]))
    ax.plot(x,y, color='b')
    ax.plot(x, uy, color='b', linestyle='dashed')
    ax.plot(x, ly, color='b', linestyle='dashed')
    ax.set_ylabel(r'predicted $\log SFR [M_{\odot} yr^{-1}]$')
    ax.set_xlabel(r'predicted $\log M_{*} [M_{\odot}]$')
    [l.set_rotation(45) for l in ax.get_xticklabels()]
    [l.set_rotation(45) for l in ax.get_yticklabels()]
    ax.set_xlim((8,12))
    ax.set_ylim((-2.5,2.0))
ax5 = P.subplot(1,5,1, aspect='auto')
triangle.hist2d(avg_mass, avg_sfr, ax=ax5, bins=50, extent=([8.0, 12.0],[-2.5,2.0]))
ax5.plot(x,y, color='b')
ax5.plot(x, uy, color='b', linestyle='dashed')
ax5.plot(x, ly, color='b', linestyle='dashed')
ax5.set_xlabel(r'$\log M_{*} [M_{\odot}]$')
ax5.set_ylabel(r'$\log SFR [M_{\odot} yr^{-1}]$')
[l.set_rotation(45) for l in ax5.get_xticklabels()]
[l.set_rotation(45) for l in ax5.get_yticklabels()]
ax5.text(8.3, 1.5, 'GZ2 galaxies')
#ax = P.subplot(1, len(past), 1)
#ax.set_ylabel(r'predicted $NUV-u$ colour')
P.tight_layout()
fig.savefig('sfr_mass_evo.pdf')

        