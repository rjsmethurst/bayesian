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
import random

font = {'family':'serif', 'size':16}
P.rc('font', **font)
P.rc('xtick', labelsize='medium')
P.rc('ytick', labelsize='medium')
P.rc('axes', labelsize='x-large')

#reason = str(raw_input('Why are you running this iteration? : '))
#
##Using PyFits to open up the Galaxy Zoo data
#file = '/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match_GZ1_k_correct_green_valley.fits'
#dat = F.open(file)
#gz2data = dat[1].data
#dat.close()
#
#col = N.zeros(14*len(gz2data)).reshape(len(gz2data), 14)
#col[:,0] = gz2data.field('MU_MR')
#col[:,1] = gz2data.field('NUV_U')
#col[:,2] = gz2data.field('t08_odd_feature_a24_merger_debiased')
#col[:,3] = gz2data.field('t06_odd_a15_no_debiased')
#col[:,4] = gz2data.field('Err_MU_MR')
#col[:,5] = gz2data.field('Err_NUV_U')
#col[:,6] = gz2data.field('z_1')
#col[:,7] = gz2data.field('zErr_1')
#col[:,8] = gz2data.field('GV_first')
#col[:,9] = gz2data.field('GV_sec')
#col[:,10] = gz2data.field('upper_GV')
#col[:,11] = gz2data.field('lower_GV')
#col[:,12] = gz2data.field('t08_odd_feature_a24_merger_count')
#col[:,13] = gz2data.field('t06_odd_a15_no_count')
#
#non_nan = N.logical_not(N.isnan(col[:,1])).astype(int)
#colours = N.compress(non_nan, col, axis=0)
#gvf = colours[colours[:,8]==1]
#colours = colours[colours[:,0] < 6.0]
#colours = colours[colours[:,1] > -3.0]
#
#mergers_count = colours[colours[:,12] > 10]
#no_mergers_count = colours[colours[:,13] > 10]
#idx = random.sample(range(len(no_mergers_count)),len(mergers_count))
#no_m_count = no_mergers_count[idx]
#merger_sample = N.append(mergers_count, no_m_count, axis=0)
#print len(merger_sample)
#
#P.figure()
#P.hist(merger_sample[:,2])
#P.show()
#P.figure()
#P.hist(merger_sample[:,3])
#P.show()
#
#
#
#age_save = '/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/mergers/age_merger_sample.npy'
#age_path = os.path.exists(age_save)
#cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
#if age_path ==False:
#    age = N.array(cosmo.age(merger_sample[:,6]))
#    N.save(age_save, age)
#else:
#    age = N.load(age_save)
#print len(age)
#
#w = [7.5, 1.5, 7.0, 1.5, 4.0, 1.5, 4.0, 1.5]
#nwalkers = 100
#nsteps= 150
#start = [7.5, 1.5, 7.5, 1.5]
#
#f = open('/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/log.txt', 'a')
#f.write('Run started at '+str(time.strftime('%H:%M'))+' on '+str(time.strftime('%d/%m/%y'))+'\n')
#f.write('Reason for running iteration: ' +reason+'\n')
#f.write('Number of walkers : '+str(nwalkers)+'\n')
#f.write('Number of steps :'+str(nsteps)+'\n')
#f.write('Starting point of walkers : '+str(start)+'\n')
#f.write('Prior assumptions, w :'+str(w)+'\n')
#f.write('Number of galaxies used was : '+str(len(age))+'\n')
#f.write('Ages of galaxies used found at : '+str(age_save)+'\n')
#f.close()
#
## w is the prior conditions on my theta parameters - the mean and standard deviation of tq and tau for the disc and smooth populations
## w = [mu_tqs, mu_taus, mu_tqd, mu_taud, sig_tqs, sig_taus, sig_tqd, sig_taud]
#start_time = time.time()
##The rest calls the emcee code and makes plots....
#samples, fig, samples_save = sample(4, nwalkers, nsteps, start, w, merger_sample[:,0], merger_sample[:,4], merger_sample[:, 1], merger_sample[:, 5], age, merger_sample[:,3], merger_sample[:,2])
#elap = (time.time() - start_time)/60
#print 'Minutes taken for '+str(len(samples)/nwalkers
#                               )+' steps and '+str(nwalkers)+' walkers', elap

samples = N.load('samples_mergers_15000_7196_16_20_08_04_14.npy')
#tqs_mcmc, taus_mcmc, tqd_mcmc, taud_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*N.percentile(samples, [16,50,84],axis=0)))
#
#f = open('/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/log.txt', 'a')
#f.write('Samples found at : '+str(samples_save)+'\n')
#f.write('Total number of positions for each parameter : '+str(len(samples))+'\n')
#f.write('t_bar from MCMC : '+str(tqs_mcmc)+'\n')
#f.write('tau_bar from MCMC : '+str(taus_mcmc)+'\n')
#f.write('t_no_bar from MCMC : '+str(tqd_mcmc)+'\n')
#f.write('tau_no_bar from MCMC : '+str(taud_mcmc)+'\n')
#f.write('Time taken to complete : '+str(elap/60)+' hours \n')
#f.write(' \n')
#f.write('------------------------------------------------------------------ \n')
#f.write(' \n')
#f.close()
#
#
#print 'tq_bar',tqs_mcmc
#print 'tau_bar',taus_mcmc
#print 'tq_no bar',tqd_mcmc
#print 'tau_no bar',taud_mcmc

fig_s = corner_plot(samples[:,0:2], labels = [r'$ t_{merger}$', r'$ \tau_{merger}$'])
fig_s.savefig('triangle_t_tau_mergers_'+str(len(samples))+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf')
fig_d = corner_plot(samples[:,2:4], labels = [r'$t_{no merger}$', r'$\tau_{no merger}$'])
fig_d.savefig('triangle_t_tau_no_mergers_'+str(len(samples))+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf')

#fig_samp = walker_plot(samples, nwalkers, nsteps)
#fig_2 = walker_steps(samples, nwalkers, nsteps)

