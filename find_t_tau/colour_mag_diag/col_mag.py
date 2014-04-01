from t_tau_func import *
import numpy as N
from matplotlib import pyplot as PY
import pylab as PY
import pyfits as F
from scipy.stats import norm, gumbel_r, gumbel_l
from scipy import linalg
import os
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
import time


font = {'family':'serif', 'size':12}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')


#Using PyFits to open up the Galaxy Zoo data
file = '/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match_GZ1_k_correct_MPA_JHU_SFR_mass.fits'
dat = F.open(file)
gz2data = dat[1].data
dat.close()
col = N.zeros(11*len(gz2data)).reshape(len(gz2data), 11)
col[:,0] = gz2data.field('MU_MR')
col[:,1] = gz2data.field('NUV_U')
col[:,2] = gz2data.field('t01_smooth_or_features_a01_smooth_debiased')
col[:,3] = gz2data.field('t01_smooth_or_features_a02_features_or_disk_debiased')
col[:,4] = gz2data.field('Err_MU_MR')
col[:,5] = gz2data.field('Err_NUV_U')
col[:,6] = gz2data.field('z_1')
col[:,7] = gz2data.field('zErr_1')
col[:,8] = gz2data.field('AVG_MASS')
col[:,9] = gz2data.field('AVG_SFR')
col[:,10] = gz2data.field('MR')
non_nan = N.logical_not(N.isnan(col[:,1])).astype(int)
colours = N.compress(non_nan, col, axis=0)
colours = colours[colours[:,8] > -1]

meanz = N.mean(colours[:,6])
cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
meanage = cosmo.age(meanz)
#age_save = '/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/colour_mag_diag/age_masses.npy'
#age_path = os.path.exists(age_save)
#cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
#if age_path ==False:
#    age = N.array(cosmo.age(colours[:,6]))
#    N.save(age_save, age)
#else:
#    age = N.load(age_save)
#print len(age)

tq = N.linspace(0.0, 13.8, 50)
tau = N.linspace(0.001, 3, 50)
tqs = N.outer(N.ones(len(tq)), tq)
taus = N.outer(tau, N.ones(len(tau)))

time = N.arange(0, 0.01, 0.003)
t = N.linspace(0, 13.7, 100)
t = N.append(time, t[1:])

ur = N.zeros((len(tau), len(tq)))
nuv = N.zeros_like(ur)
rmag = N.zeros_like(ur)
mass = N.zeros_like(ur)
k = 0
for n in range(len(tq)):
    for m in range(len(tau)):
        ur[m,n], nuv[m,n], rmag[m,n], mass[m,n] = predict_c_one([tq[n], tau[m]], meanage)
        k+=1
        print k

N.save('ur.npy', ur)
N.save('nuv.npy', nuv)
N.save('rmag.npy', rmag)
N.save('mass.npy', mass)

#ur = N.load('ur.npy')
#nuv = N.load('nuv.npy')
#rmag = N.load('rmag.npy')
#mass = N.load('mass.npy')


P.figure()
ax=P.subplot(121)
triangle.hist2d(colours[:,8], colours[:,0], ax=ax, extent=[(8,12),(0.5,3.5)], bins=50, plot_datapoints=False)
[l.set_rotation(45) for l in ax.get_xticklabels()]
[l.set_rotation(45) for l in ax.get_yticklabels()]
ax.set_xlabel(r'$M_r$')
ax.set_ylabel(r'$u-r$')
ax1 = P.subplot(122)
ax1.scatter(N.log10(mass), ur, c=tqs, s=taus*10)
[l.set_rotation(45) for l in ax1.get_xticklabels()]
[l.set_rotation(45) for l in ax1.get_yticklabels()]
ax1.set_xlabel(r'$M_r$')
ax1.set_ylabel(r'$u-r$')
P.tight_layout()
P.savefig('col_mag_diagram.pdf')
P.show()

