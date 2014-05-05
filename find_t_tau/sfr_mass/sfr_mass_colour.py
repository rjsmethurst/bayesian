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
#
#
#ur = N.load('ur_new_sfhs_avg_age.npy')
#nuv = N.load('nuvu_new_sfhs_avg_age.npy')

#N.save('sfr_new_sfhs_avg_age.npy', sfr)
#N.save('mass_new_sfhs_avg_age.npy', mass)
#
#N.save('full_sfr.npy', full_sfr)
#N.save('full_mass.npy', full_mass)

sfr = N.load('sfr_new_sfhs_avg_age.npy')
mass = N.load('mass_new_sfhs_avg_age.npy')

#P.figure(figsize=(30,5))
#for n in range(len(age)):
#    ax = P.subplot(1, len(age)+1, n+2)
#    triangle.hist2d(ur[:,n], nuv[:,n], ax=ax,bins=50, extent=([0.0, 3.0],[-1.0,5.0]))
#    [l.set_rotation(45) for l in ax.get_xticklabels()]
#    [l.set_rotation(45) for l in ax.get_yticklabels()]
#    ax.set_ylabel(r'predicted $NUV-u$')
#    ax.set_xlabel(r'predicted $u-r$')
#    ax.text(0.25, 4.5, 't = '+str(age[n])+' Gyr')
#    print n
#ax5 = P.subplot(1,6,1)
#triangle.hist2d(u_r, nuv_u, ax=ax5, bins=50, extent=([0.0, 3.0],[-1.0,5.0]))
#ax5.set_xlabel(r'$u-r$')
#ax5.set_ylabel(r'$NUV-u$')
#[l.set_rotation(45) for l in ax5.get_xticklabels()]
#[l.set_rotation(45) for l in ax5.get_yticklabels()]
#ax5.text(0.25, 4.5, 'GZ2 galaxies')
#P.tight_layout()
#P.savefig('pred_colour.pdf')
#P.show()

def peng_sfr(m,t):
    return (2.5*((m/10**10)**(-0.1))*((t/3.5)**(-2.2)))*(ms/1E9)

ms = N.array((1E8, 1E12, 51))
#m = (1.0--1)/(11.0-9)
#y = m*x + (-1-(9*m))
#uy = y+0.3
#ly = y-0.3
sfr_138 = N.log10(peng_sfr(ms, 13.8))
sfr_12 = N.log10(peng_sfr(ms, 12.0))
sfr_9 = N.log10(peng_sfr(ms, 10.0))
sfr_6 = N.log10(peng_sfr(ms, 8.0))
sfr_3 = N.log10(peng_sfr(ms, 6.0))
m = N.log10(ms)
             
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

fig = P.figure(figsize=(13,9))

ax1 = P.subplot(2,3,1, aspect='auto')
triangle.hist2d(avg_mass, avg_sfr, ax=ax1, bins=50, extent=([8.0, 12.0],[-2.5,2.0]))
ax1.plot(m,sfr_138, color='b')
ax1.plot(m, sfr_138+0.3, color='b', linestyle='dashed')
ax1.plot(m, sfr_138-0.3, color='b', linestyle='dashed')
ax1.tick_params(labelbottom='off')
ax1.set_yticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])
ax1.set_yticklabels([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])
[l.set_rotation(45) for l in ax1.get_xticklabels()]
[l.set_rotation(45) for l in ax1.get_yticklabels()]
ax1.set_ylabel(r'$\log SFR [M_{\odot} yr^{-1}]$')
ax1.set_ylim(-2.5, 2.0)
ax1.text(8.3, 1.5, 'GZ2 galaxies')

ax2 = P.subplot(2,3,2, aspect='auto')
triangle.hist2d(N.log10(mass[:,0]), log_sfr[:,0], ax=ax2, bins=50, extent=([8.0, 12.0],[-2.5, N.max(log_sfr[:,0])]))
ax2.plot(m,sfr_138, color='b')
ax2.plot(m, sfr_138+0.3, color='b', linestyle='dashed')
ax2.plot(m, sfr_138-0.3, color='b', linestyle='dashed')
ax2.tick_params(labelbottom='off')
ax2.tick_params(labelleft='off')
ax2.set_ylim(-2.5, 2.0)
ax2.text(8.3, 1.5, 't = '+str(age[0])+' Gyr')


ax3 = P.subplot(2,3,3, aspect='auto')
triangle.hist2d(N.log10(mass[:,1]), log_sfr[:,1], ax=ax3, bins=50, extent=([8.0, 12.0],[-2.5, N.max(log_sfr[:,1])]))
ax3.plot(m,sfr_9, color='b')
ax3.plot(m, sfr_9+0.3, color='b', linestyle='dashed')
ax3.plot(m, sfr_9-0.3, color='b', linestyle='dashed')
ax3.tick_params(labelbottom='off')
ax3.tick_params(labelleft='off')
ax3.text(8.3, 1.5, 't = '+str(age[1])+' Gyr')
ax3.set_ylim(-2.5, 2.0)

ax4 = P.subplot(2,3,4, aspect='auto')
triangle.hist2d(N.log10(mass[:,2]), log_sfr[:,2], ax=ax4, bins=50, extent=([8.0, 12.0],[-2.5, N.max(log_sfr[:,2])]))
ax4.plot(m,sfr_12, color='b')
ax4.plot(m, sfr_12+0.3, color='b', linestyle='dashed')
ax4.plot(m, sfr_12-0.3, color='b', linestyle='dashed')
ax4.set_yticks([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
ax4.set_yticklabels([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
ax4.set_xticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5])
ax4.set_xticklabels([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5])
[l.set_rotation(45) for l in ax4.get_xticklabels()]
[l.set_rotation(45) for l in ax4.get_yticklabels()]
ax4.set_ylabel(r'$\log SFR [M_{\odot} yr^{-1}]$')
ax4.text(8.3, 1.5, 't = '+str(age[2])+' Gyr')
ax4.set_ylim(-2.5, 2.0)


ax5 = P.subplot(2,3,5, aspect='auto')
triangle.hist2d(N.log10(mass[:,3]), log_sfr[:,3], ax=ax5, bins=50, extent=([8.0, 12.0],[-2.5, N.max(log_sfr[:,3])]))
ax5.plot(m,sfr_6, color='b')
ax5.plot(m, sfr_6+0.3, color='b', linestyle='dashed')
ax5.plot(m, sfr_6-0.3, color='b', linestyle='dashed')
ax5.tick_params(labelleft='off')
ax5.set_xticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5])
ax5.set_xticklabels([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5])
[l.set_rotation(45) for l in ax5.get_xticklabels()]
ax5.set_xlabel(r'$\log M_{*} [M_{\odot}]$')
ax5.text(8.3, 1.5, 't = '+str(age[3])+' Gyr')
ax5.set_ylim(-2.5, 2.0)

ax6 = P.subplot(2,3,6, aspect='auto')
triangle.hist2d(N.log10(mass[:,4]), log_sfr[:,4], ax=ax6, bins=50, extent=([8.0, 12.0],[-2.5, N.max(log_sfr[:,4])]))
ax6.plot(m,sfr_3, color='b')
ax6.plot(m, sfr_3+0.3, color='b', linestyle='dashed')
ax6.plot(m, sfr_3-0.3, color='b', linestyle='dashed')
ax6.set_xticks([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
ax6.set_xticklabels([8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0])
[l.set_rotation(45) for l in ax6.get_xticklabels()]
ax6.tick_params(labelleft='off')
ax6.text(8.3, 1.5, 't = '+str(age[4])+' Gyr')
ax6.set_ylim(-2.5, 2.0)

fig.tight_layout()
fig.subplots_adjust(wspace=0.0, hspace=0.0)
fig.savefig('sfr_mass_evo_second.pdf')

        