import numpy as N
from matplotlib import pyplot as PY
import pylab as P
import pyfits as F
import os
from astropy.cosmology import FlatLambdaCDM
import time
import triangle



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
col = N.zeros(15*len(gz2data)).reshape(len(gz2data), 15)
col[:,0] = gz2data.field('MU_MR')
col[:,1] = gz2data.field('NUV_U')
col[:,2] = gz2data.field('t01_smooth_or_features_a01_smooth_debiased')
col[:,3] = gz2data.field('t01_smooth_or_features_a02_features_or_disk_debiased')
col[:,4] = gz2data.field('Err_MU_MR')
col[:,5] = gz2data.field('Err_NUV_U')
col[:,6] = gz2data.field('z_1')
col[:,7] = gz2data.field('zErr_1')
col[:,8] = gz2data.field('upper_GV')
col[:,9] = gz2data.field('lower_GV')
col[:,10] = gz2data.field('MR')
col[:,11] = gz2data.field('t03_bar_a06_bar_count')
col[:,12] = gz2data.field('t03_bar_a07_no_bar_count')
col[:,13] = gz2data.field('t03_bar_a06_bar_debiased')
col[:,14] = gz2data.field('t03_bar_a07_no_bar_debiased')

non_nan = N.logical_not(N.isnan(col[:,1])).astype(int)
colours = N.compress(non_nan, col, axis=0)
colours = colours[colours[:,8] > -1]

bars_count = colours[colours[:,12] > 10]
no_bars_count = colours[colours[:,13] > 10]
bar_sample = N.append(bars_count, no_bars_count, axis=0)
print len(bar_sample)
bars = bar_sample[bar_sample[:,13] > 0.5]
nobars = bar_sample[bar_sample[:,14] > 0.5]

#gv1 = colours[colours[:,0] < colours[:,8]]
#gv = gv1[gv1[:,0] > gv1[:,9]]
#rs = red_s = colours[colours[:,0] > colours[:,8]]
#bc = blue_c = colours[colours[:,0] < colours[:,9]]
#
#gv_disc = gv[gv[:,3] > 0.5]
#gv_smooth = gv[gv[:,2] > 0.5]
#gv_late = gv[gv[:,3] >= 0.8]
#gv_early = gv[gv[:,2] >=0.8]
#
#bc_disc = bc[bc[:,3] > 0.5]
#bc_smooth = bc[bc[:,2] > 0.5]
#bc_late = bc[bc[:,3] >= 0.8]
#bc_early = bc[bc[:,2] >=0.8]
#
#rs_disc = rs[rs[:,3] > 0.5]
#rs_smooth = rs[rs[:,2] > 0.5]
#rs_late = rs[rs[:,3] >= 0.8]
#rs_early = rs[rs[:,2] >=0.8]



P.figure(figsize=(13,5))
ax1=P.subplot(1, 3, 1)
triangle.hist2d(colours[:,0], colours[:,1], ax=ax1, extent=[(0.75,3.25), (-1, 5)], bins=75, plot_datapoints=True)
[l.set_rotation(45) for l in ax1.get_xticklabels()]
ax1.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax1.set_xticklabels([1.0, 1.5, 2.0, 2.5, 3.0])
ax1.set_ylabel(r'$NUV - u$')
ax1.text(0.85, 4.5, 'All GZ2 galaxies', color='k')
[l.set_rotation(45) for l in ax1.get_yticklabels()]

ax2 = P.subplot(1,3,2, sharey=ax1, sharex=ax1)
triangle.hist2d(bars[:,0], bars[:,1], ax=ax2, extent=[(0.75,3.25), (-1, 5)], bins=40, plot_datapoints=True)
[l.set_rotation(45) for l in ax2.get_xticklabels()]
ax2.tick_params(axis='y', labelleft='off')
ax2.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax2.set_xticklabels([1.0, 1.5, 2.0, 2.5, 3.0])
ax2.set_xlabel(r'$u-r$')
ax2.text(0.85, 4.5, 'Barred disc galaxies', color='k')

ax3 = P.subplot(1,3,3, sharey=ax1, sharex = ax1)
triangle.hist2d(nobars[:,0], nobars[:,1], ax=ax3, extent=[(0.75,3.25), (-1, 5)], bins=50, plot_datapoints=True)
[l.set_rotation(45) for l in ax3.get_xticklabels()]
ax3.tick_params(axis='y', labelleft='off')
ax3.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax3.set_xticklabels([1.0, 1.5, 2.0, 2.5, 3.0])
ax3.text(0.85, 4.5, 'Unbarred disc galaxies', color='k')
P.tight_layout()
P.subplots_adjust(wspace=0.0)
P.savefig('bar_no_bars_col_col.pdf')
P.show()



#P.figure(figsize=(13,9))
#ax1=P.subplot(2, 3, 1)
#triangle.hist2d(gv[:,0], gv[:,1], ax=ax1, extent=[(0.75,3.25), (-1, 5)], bins=75, plot_contours=True)
#[l.set_rotation(45) for l in ax1.get_xticklabels()]
#ax1.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
#ax1.set_xticklabels([1.0, 1.5, 2.0, 2.5, 3.0])
#ax1.set_yticks([-1, 0, 1, 2, 3, 4, 5])
#ax1.set_yticklabels([-1, 0, 1, 2, 3, 4, 5])
#ax1.set_ylabel(r'$NUV - u$')
#ax1.text(0.85, 4.5, 'All Green Valley', color='Black')
#[l.set_rotation(45) for l in ax1.get_yticklabels()]
#
#ax2 = P.subplot(2,3,2, sharey=ax1, sharex=ax1)
#triangle.hist2d(gv_smooth[:,0], gv_smooth[:,1], ax=ax2, extent=[(0.75,3.25), (-1, 5)], bins=50, plot_contours=True, color = 'Red', linewidths=1.5)
#ax2.tick_params(axis='y', labelleft='off')
#ax2.tick_params(axis='x', labelbottom='off')
#ax2.text(0.85, 4.5, 'Smooth-like', color='Red')
#
#ax3 = P.subplot(2,3,3, sharey=ax1, sharex = ax1)
#triangle.hist2d(gv_early[:,0], gv_early[:,1], ax=ax3, extent=[(0.75,3.25), (-1, 5)], bins=50, plot_contours=True, color='Red', linewidths=1.5)
#ax3.tick_params(axis = 'x', labelbottom='off')
#ax3.tick_params(axis='y', labelleft='off')
#ax3.text(0.85, 4.5, 'Clean Early-tpyes', color='Red')
#
#ax4 = P.subplot(2,3,5, sharey=ax1, sharex=ax1)
#triangle.hist2d(gv_disc[:,0], gv_disc[:,1], ax=ax4, extent=[(0.75,3.25), (-1, 5)], bins=50, plot_contours=True, color = 'Blue', linewidths=1.5)
#[l.set_rotation(45) for l in ax4.get_xticklabels()]
#[l.set_rotation(45) for l in ax4.get_yticklabels()]
#ax4.set_ylabel(r'$NUV-u$')
#ax4.set_yticks([-1, 0, 1, 2, 3, 4])
#ax4.set_yticklabels([-1, 0, 1, 2, 3, 4])
#ax4.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
#ax4.set_xticklabels([1.0, 1.5, 2.0, 2.5, 3.0])
#ax4.set_xlabel(r'$u-r$')
#ax4.text(0.85, 4.5, 'Disc-like', color='Blue')
#
#
#ax5 = P.subplot(2,3,6, sharey=ax1, sharex=ax1)
#triangle.hist2d(gv_late[:,0], gv_late[:,1], ax=ax5, extent=[(0.75,3.25), (-1, 5)], bins=50, plot_contours=True, color='Blue', linewidths=1.5)
#[l.set_rotation(45) for l in ax5.get_xticklabels()]
#ax5.tick_params(axis='y', labelleft='off')
#ax5.set_xticks([1.0, 1.5, 2.0, 2.5, 3.0])
#ax5.set_xticklabels([1.0, 1.5, 2.0, 2.5, 3.0])
#ax5.text(0.85, 4.5, 'Clean Late-types', color='Blue')
#ax5.set_xlabel(r'$u-r$')
#
#P.tight_layout()
#P.subplots_adjust(wspace=0.0, hspace=0.0)
#P.savefig('gv_clean_col_col.pdf')
#P.show()

