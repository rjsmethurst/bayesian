import t_tau_function
import numpy as N
from matplotlib import pyplot as PY
import pylab as PY
from hist_function import *
import pyfits as F
from scipy.stats import norm, gumbel_r, gumbel_l
from scipy import linalg
import os
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM

font = {'family':'serif', 'size':12}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')

file = '/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match.fits'
dat = F.open(file)
gz2data = dat[1].data

col = N.zeros(19*len(gz2data)).reshape(len(gz2data), 19)
col[:,0] = gz2data.field('petroMag_u') - gz2data.field('petroMag_r')
col[:,1] = gz2data.field('NUV') - gz2data.field('petroMag_u')
col[:,2] = gz2data.field('t01_smooth_or_features_a01_smooth_flag')
col[:,3] = gz2data.field('t01_smooth_or_features_a02_features_or_disk_flag')
#col[:,4] = gz2data.field('t03_bar_a06_bar_flag')
#col[:,5] = gz2data.field('t03_bar_a06_bar_debiased')
#col[:,6] = gz2data.field('t04_spiral_a08_spiral_flag')
#col[:,7] = gz2data.field('t11_arms_number_a31_1_debiased')
#col[:,8] = gz2data.field('t11_arms_number_a32_2_debiased')
#col[:,9] = gz2data.field('t11_arms_number_a33_3_debiased')
#col[:,10] = gz2data.field('t11_arms_number_a34_4_debiased')
#col[:,11] = gz2data.field('t11_arms_number_a36_more_than_4_debiased')
#col[:,12] = gz2data.field('t11_arms_number_a37_cant_tell_debiased')
#col[:,13] = gz2data.field('t10_arms_winding_a28_tight_debiased')
#col[:,14] = gz2data.field('t10_arms_winding_a29_medium_debiased')
#col[:,15] = gz2data.field('t10_arms_winding_a30_loose_flag')
col[:,16] = ((gz2data.field('petroMagErr_u')**2 + gz2data.field('petroMagErr_r')**2)**0.5)
col[:,17] = ((gz2data.field('petroMagErr_u')**2+ gz2data.field('e_NUV')**2)**0.5)
col[:,18] = gz2data.field('z')
 

non_nan = N.logical_not(N.isnan(col[:,1])).astype(int)
colours = N.compress(non_nan, col, axis=0)

disc = colours[colours[:,3]==1]
smooth = colours[colours[:,2]==1]
print len(smooth)
#inter = colours[colours[:,2]!=1]
#intermediate = inter[inter[:,3]!=1]
#spiral = colours[colours[:,6]==1]
#bars = colours[colours[:,5]>=0.5]

#Define 'made-up' test galaxy properties 
#z_gal = N.array([0.15, 0.21])
#u_r_gal = N.array([1.6, 2.1])
#u_r_gal_sigma = N.array([0.01, 0.015])
#nuv_u_gal = N.array([0.5, 1.2])
#nuv_u_gal_sigma = N.array([0.02, 0.015])
u_r_gal = smooth[100:101,0]
print N.shape(u_r_gal)
print 'u_r_gal', u_r_gal
u_r_gal_sigma = smooth[100:101,16]
print 'u_r_gal_sigma', u_r_gal_sigma
nuv_u_gal = smooth[100:101,1]
nuv_u_gal_sigma = smooth[100:101,17]

age_save = '/Users/becky/Projects/Green-Valley-Project/bayesian/age_smooth.npy'
age_path = os.path.exists(age_save)
cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
if age_path == True:
    age = N.array(cosmo.age(smooth[100:101,18]))
    N.save(age_save, age)
else:
    age = N.load(age_save)
#age = cosmo.age(z_gal)
print age
print N.shape(age)
#
#col_gal = N.array([u_r_gal, nuv_u_gal])
#sigma_gal = N.array([u_r_gal_sigma, nuv_u_gal_sigma])

#Define model properties
#time = (N.arange(28.0)/2.0) # In Gyrs
#time = N.linspace(0.0, cosmo.age(0.0), 50) # In Gyrs
#tim = N.arange(0, 0.01, 0.003)
#times = N.arange(0, 13.7, 0.01)
#t = N.append(tim, times[1:])
#tau = N.array([0.001, 0.25, 1.0, 2.5])
#tau = N.linspace(0.001, 2.5, 50)
#taus = N.outer(tau, N.ones(len(tau)))
##tq = N.array([8,9,10])
#tq = N.linspace(5.0, 11.0, 50)
#tqs = N.outer(N.ones(len(tq)), tq)
#ext = (N.min(tq), N.max(tq), N.min(tau), N.max(tau))
#
## Predict colours for test galaxy for each value of tau and tq either by reading an already saved file or creating it using the t_tau_function module
#nuv_u_save = '/Users/becky/Projects/Green-Valley-Project/bayesian/nuv_u_pred_'+str(len(tau))+'_'+str(len(tq))+'_'+str(len(time))+'.npy'
#u_r_save = '/Users/becky/Projects/Green-Valley-Project/bayesian/u_r_pred_'+str(len(tau))+'_'+str(len(tq))+'_'+str(len(time))+'.npy'
#nuv_u_path = os.path.exists(nuv_u_save)
#u_r_path = os.path.exists(u_r_save)
#if nuv_u_path == False or u_r_path == False:
#    nuv_u_pred, u_r_pred = t_tau_function.predict_colour_numpy(tq, tau, time)
#    N.save(nuv_u_save, nuv_u_pred)
#    N.save(u_r_save, u_r_pred)
#else:
#    nuv_u_pred = N.load(nuv_u_save)
#    u_r_pred = N.load(u_r_save)
#
#print u_r_pred[0:2,47:49,0:2]
#
##Interp between the time values to find the predicted colours for the given age(z) of the galaxy at all t and tau
#nuv_u_ages = t_tau_function.age_colours(nuv_u_pred, tau, time, tq, age)
#u_r_ages = t_tau_function.age_colours(u_r_pred, tau, time, tq, age)
#print N.shape(nuv_u_ages)
#print N.shape(u_r_ages)
#print 'u_r_ages',u_r_ages[0:2,0,0:2]

#Calculate log_p_d_giv_theta which is an array of dim (len(tq), len(tau)). The predicted colours are transposed due to the shape of the array when spliced
#log_p_d_giv_theta_u_r = t_tau_function.log_p_d_given_theta(u_r_gal, u_r_gal_sigma, N.array(N.split(u_r_ages.T, len(age), axis=1)))
#log_p_d_giv_theta_nuv_u = t_tau_function.log_p_d_given_theta(nuv_u_gal, nuv_u_gal_sigma, N.array(N.split(nuv_u_ages.T, len(age), axis=1)))
#
#print N.shape(log_p_d_giv_theta_u_r)
#print 'log(P(d|theta))_u_r', log_p_d_giv_theta_u_r[0:2,0:2]
#
##Define w with mu_tq = w[0],  mu_tau = w[1], sig_tq =  w[2], sig_tau = w[3]
w = [7.0, 0.5, 2.35, 0.7]
#log_prob = N.zeros(len(tau)*len(tq)).reshape(len(tq), len(tau))
#for n in range(len(tau)):
#    for m in range(len(tq)):
#        theta = N.array([tq[m], tau[n]]).reshape(2,1)
#        log_prob[m,n] = t_tau_function.log_p_theta(theta, w)
#
#print 'sum log prob', N.sum(N.exp(log_prob))
#
#dist_tq = norm(loc=w[0], scale=w[2])
#dist_tau = norm(loc=w[1], scale = w[3])
##
#pdf_tau = dist_tau.pdf(tau)
#pdf_tq = dist_tq.pdf(tq)
#
#sum = log_p_d_giv_theta_u_r + log_prob
#print 'sum', sum[0:2,0:2]
#p_d_theta_p_theta_u_r = N.exp(sum)
#print 'P(d|theta)P(theta)', p_d_theta_p_theta_u_r[0:2,0:2]
#norm_factor_u_r = N.sum(p_d_theta_p_theta_u_r)
#print 'norm', norm_factor_u_r
#p_theta_given_data_u_r = p_d_theta_p_theta_u_r/norm_factor_u_r
#print 'P(theta|d)', p_theta_given_data_u_r[0:2,0:2]
#print N.shape(p_theta_given_data_u_r)
#print N.sum(p_theta_given_data_u_r)
#
#p_d_theta_p_theta_nuv_u = N.exp(log_p_d_giv_theta_nuv_u + log_prob)
#norm_factor_nuv_u = N.sum(p_d_theta_p_theta_nuv_u)
#p_theta_given_data_nuv_u = p_d_theta_p_theta_nuv_u/norm_factor_nuv_u
#print N.shape(p_theta_given_data_nuv_u)
#print p_theta_given_data_nuv_u
#print N.sum(p_theta_given_data_nuv_u)

#lev = N.array([0.1, 0.3, 0.5, 0.7, 0.9])
#levels = lev*N.max(p_theta_given_data_nuv_u)
#print N.max(p_theta_given_data_nuv_u)
#print levels
#
#
#ax1 = PY.subplot(2,2,2)
#cont = ax1.contour(log_prob, extent=ext, cmap=PY.cm.binary)
#ax1.set_ylabel(r'$\tau  (Gyr)$')
#ax1.set_xlabel(r'$t_{quench}  (Gyr)$')
#ax2 = PY.subplot(2,2,4)
#ptq = P.plot(tq, pdf_tq, color='k')
#ptq = PY.gca()
#ptq.invert_yaxis()
#ax2.tick_params(axis='x', labeltop='off', labelbottom='off')
#ax2.set_ylabel('Probability Density')
#ax3 = PY.subplot(2,2,1)
#ptau = PY.plot(pdf_tau, tau, color='k')
#ptau = PY.gca()
#ptau.invert_xaxis()
#ax3.tick_params(axis='y', labelright='off', labelleft='off')
#ax3.set_xlabel('Probability Density')
#PY.subplots_adjust(wspace = 0.25, hspace =0.25)
##PY.title('$P(\theta|w)$')
#PY.show()
#
#P.figure()
#ax1 = P.contour(log_p_d_giv_theta_nuv_u, origin = 'upper', extent=ext, cmap=P.cm.binary)
#P.ylabel(r'$\tau  (Gyr)$')
#P.xlabel(r'$t_{quench}  (Gyr)$')
#P.title(r'$P(d(u-r)|\theta)$')
#P.colorbar()
#P.show()
#
#P.figure()
#ax1 = P.imshow(log_p_d_giv_theta_nuv_u, extent=ext, aspect='auto', cmap=P.cm.binary)
#P.ylabel(r'$\tau  (Gyr)$')
#P.xlabel(r'$t_{quench}  (Gyr)$')
#P.title(r'$P(d(u-r)|\theta)$')
#P.colorbar()
#P.show()
#
#P.figure()
#ax1 = P.contour(p_theta_given_data_nuv_u, origin='upper', extent=ext, cmap=P.cm.binary)
#P.ylabel(r'$\tau  (Gyr)$')
#P.xlabel(r'$t_{quench}  (Gyr)$')
#P.title(r'$P(\theta|d(u-r))$')
#P.colorbar()
#P.show()
#
#P.figure()
#ax1 = P.imshow(p_theta_given_data_nuv_u, extent=ext, aspect='auto', cmap=P.cm.binary)
#P.ylabel(r'$\tau  (Gyr)$')
#P.xlabel(r'$t_{quench}  (Gyr)$')
#P.title(r'$P(\theta|d(u-r))$')
#P.colorbar()
#P.show()

               
samples, fig = t_tau_function.sample(2, 64, u_r_gal, u_r_gal_sigma, age, w)
tq_mcmc, tau_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*N.percentile(samples, [16,50,84],axis=0)))
print tq_mcmc
print tau_mcmc
