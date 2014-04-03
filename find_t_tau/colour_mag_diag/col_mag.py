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


font = {'family':'serif', 'size':14}
P.rc('font', **font)
P.rc('xtick', labelsize='medium')
P.rc('ytick', labelsize='medium')
P.rc('axes', labelsize='large')


#Using PyFits to open up the Galaxy Zoo data
file = '/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match_GZ1_k_correct_green_valley.fits'
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
col[:,10] = gz2data.field('MR')
non_nan = N.logical_not(N.isnan(col[:,1])).astype(int)
colours = N.compress(non_nan, col, axis=0)
colours = colours[colours[:,8] > -1]

Mr = N.linspace(N.min(colours[:,10]), N.max(colours[:,10]), 200)
C_dash = 2.06 - 0.244*N.tanh((Mr + 20.07)/1.09)
upper = C_dash + 0.128
lower = C_dash - 0.128


P.figure()
ax=P.subplot(111)
triangle.hist2d(colours[:,10], colours[:,0], ax=ax, extent=[(-24, -18),(0.5,3.5)], bins=75, plot_datapoints=False)
ax.plot(Mr, C_dash, 'g', ls='dashed')
ax.plot(Mr, upper, 'g', lw = 2)
ax.plot(Mr, lower, 'g', lw=2)
[l.set_rotation(45) for l in ax.get_xticklabels()]
[l.set_rotation(45) for l in ax.get_yticklabels()]
ax.set_xlim((-18, -24))
ax.set_xlabel(r'$M_r$')
ax.set_ylabel(r'$u-r$')
P.tight_layout()
P.savefig('col_mag_with GV.pdf')
P.show()

