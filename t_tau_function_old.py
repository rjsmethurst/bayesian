"""This function will return a predicted color value for given tau and t_quench (starting with scalars and working up to numpy arrays) calculated at the age of the galaxy (inferred from it's redshift assuming all the galaxies formed at t=0, as in the models)."""

##############################################################################
#
# BC03 models are obtained from www.iap.fr/~charlot/bc03
# The files from the various bc03.*.tar.gz archives are organised into four
# subdirectories which will be called in this code:
#
# ./bc03/doc
# ./bc03/models
# ./bc03/src
# ./bc03/templates
#
# The BC03 models are defined in a document provided in the download and
# can be found in ./bc03/doc/bc03.ps
#
### INPUTS
#
# tau
# t_quench
# redshift of galaxy - 1D array (age calculated from this using astropy cosmo calculator)
# BC03 model files
#
### OUPUTS
#
# Predicted colour for age of given galaxy and tau and tquench values provided
#
##############################################################################

import numpy as N
import scipy as S
import pylab as P
import pyfits as F
import idlsave
from scipy.integrate import simps
import scipy.spatial.distance as SD
from hist_function import *
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
import assign_fluxes
import pyfits as F
from scipy.stats import norm, gumbel_r, gumbel_l
from scipy import linalg
from itertools import product
from scipy.interpolate import griddata
import emcee
import triangle

cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)


def predict_colour_for_loop(tq, tau, z):
    cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
    age = cosmo.age(z)
    print age
    #time = N.arange(0, 0.01, 0.003)
    #t = N.arange(0, 13.7, 0.01)
    #time = N.append(time, t[1:])
    time = N.arange(28.0)/2.0
    dir ='/Users/becky/Projects/Green-Valley-Project/bc03/models/Padova1994/chabrier/ASCII/'
    model = 'extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
    data = N.loadtxt(dir+model)
    sfr = N.zeros(len(time)*len(tq)*len(tau)).reshape(len(time), len(tq), len(tau))
    nuv_u = N.zeros_like(sfr)
    u_r = N.zeros_like(sfr)
    nuv_u_age = N.zeros(len(age)*len(tq)*len(tau)).reshape(len(age),len(tq), len(tau))
    u_r_age = N.zeros_like(nuv_u_age)
    for m in range(len(tq)):
        for n in range(len(tau)):
            sfr[:,m,n] = expsfh(tau[n], tq[m], time)
            total_flux = assign_fluxes.assign_total_flux(data[0,1:], data[1:,0], data[1:,1:], time*1E9, sfr[:,m,n])
            nuv_u[:,m,n], u_r[:,m,n] = get_colours(total_flux, data)
            nuv_u_age[:,m,n] = N.interp(age, time, nuv_u[:,m,n])
            u_r_age[:,m,n] = N.interp(age, time, u_r[:,m,n])
    return nuv_u_age, u_r_age

def predict_colour_numpy(tq, tau, t):
    # Time, tq and tau are in units of Gyrs 
    time = t.reshape(len(t),1)
    tq = tq.reshape(len(tq),1)
    tau = tau.reshape(len(tau),1)
    sfr = exp_sfh(tau, tq, time)
    a = N.split(sfr,len(tq), axis=2)
    dir ='/Users/becky/Projects/Green-Valley-Project/bc03/models/Padova1994/chabrier/ASCII/'
    model = 'extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
    data = N.loadtxt(dir+model)
    nuv_u = N.zeros_like(sfr)
    u_r = N.zeros_like(sfr)
    for m in range(len(tq)):
        for n in range(len(tau)):
            total_flux = assign_fluxes.assign_total_flux_numpy(data[0,1:], data[1:,0], data[1:,1:], t*1E9, a[m][n])
            nuv_u[n,:,m], u_r[n,:,m] = get_colours_numpy(t*1E9, total_flux, data)
    return nuv_u, u_r


def exp_sfh(tau, tq, time):
    top = SD.cdist(time, tq, lambda u,v: u-v if (u>v) else 0.0)
    frac = top/(N.split(tau, len(tau)))
    sfr = 100*(N.exp(-frac))
    return sfr

def expsfh(tau, tq, time):
    a = time.searchsorted(tq)
    sfr = N.ones(len(time))*100
    sfr[a:] = 100*N.exp(-(time[a:]-tq)/tau)
    return sfr

def age_colours(colour, tau, time, tq, ages):
    ### function to take the 3d array of colours with dimensions (len(tau),len(time), len(tq)) and interpolate for the ages of the observed galaxy to give a 3d array of dimensions (len(tau), len(ages), len(tq))
    ### this is done by flattening the colours array to calculate at the given ages and then reshaping the array at the end
    print N.shape(colour)
    print N.shape(tau)
    print N.shape(time)
    print N.shape(tq)
    print N.shape(ages)
    points = N.array(list(product(tau, time, tq))) # all combinations of tau, time, tq in colours array
    print points
    print len(points)
    new_points = N.array(list(product(tau, ages, tq))) # all combinations of tau, age, tq in new array
    print new_points
    print len(new_points)
    _, tau_index = N.unique(tau, return_inverse = True) # indices of tau array
    _, time_index = N.unique(time, return_inverse = True) # indices of time array
    _, tq_index = N.unique(tq, return_inverse = True) # indices of tq array
    point_index = N.array(list(product(tau_index, time_index, tq_index))) # all combinations of indices to match combinations of tau, time and tq in points
    print N.shape(point_index)
    values = colour[point_index[:,0], point_index[:,1], point_index[:,2]] #Flatten the data array into 1D array
    print len(values)
    new_values = griddata(points, values, new_points) # set up grid for values at given points and interpolate values for new_points
    print len(new_values)
    age_colours = new_values.reshape(len(tau), len(ages), len(tq)) # reshape from 1D into 3D as required
    return age_colours
    
def log_p_d_given_theta_old(col_gal, col_sigma_gal, col_pred):
    ### col_pred should be a 2D array of dim len(tq) x len(tau) and contain a predicted colour at all of these combinations of parameters for one colour (either nuv-u or u-r)
    ### col_gal and col_sigma_gal should be an array of galaxy colours (either nuv-u or u-r)
    print N.shape(col_pred)
    print N.shape(col_gal)
    col_gal = col_gal.reshape(len(col_gal), 1)
    col_sigma_gal = col_sigma_gal.reshape(len(col_sigma_gal), 1)
    col_gals = N.split(col_gal, len(col_gal))
    col_sigma_gals = N.array(N.split(col_sigma_gal, len(col_sigma_gal)))
    print N.shape(col_gals)
    print N.shape(col_sigma_gals)
    chisq_array = ((col_gals - col_pred[:,:,0,:])**2)/(col_sigma_gals**2)
    print col_gals - col_pred[:,:,0,:]
    print ((col_gals - col_pred[:,:,0,:])**2)
    print 'chisq', chisq_array[0,:,:]
    print N.shape(chisq_array)
    #chisq = N.sum(chisq_array)
    n = (2*(N.pi)*(col_sigma_gal**2))**0.5
    norm = N.split(n, len(col_sigma_gal))
    print norm
    log_p_d_given_theta = - N.log(norm) - (chisq_array/2.0)
    print 'log p for one', log_p_d_given_theta[0,:,:]
    total_log_p_d_given_theta = N.sum(log_p_d_given_theta, axis=0)
    print total_log_p_d_given_theta
    print N.sum(total_log_p_d_given_theta)
    return total_log_p_d_given_theta


def log_p_d_given_theta_old2(col_gal, col_sigma_gal, col_pred):
    ### col_pred should be a 2D array of dim len(tq) x len(tau) and contain a predicted colour at all of these combinations of parameters for one colour (either nuv-u or u-r)
    ### col_gal and col_sigma_gal should be an array of galaxy colours (either nuv-u or u-r)
    col_gal = col_gal.reshape(len(col_gal), 1)
    col_sigma_gal = col_sigma_gal.reshape(len(col_sigma_gal), 1)
    col_gals = N.split(col_gal, len(col_gal))
    col_sigma_gals = N.array(N.split(col_sigma_gal, len(col_sigma_gal)))
    chisq_array = ((col_gals - col_pred[:,:,0,:])**2)/(col_sigma_gals**2)
    print col_gals - col_pred[:,:,0,:]
    print ((col_gals - col_pred[:,:,0,:])**2)
    print 'chisq', chisq_array[0,:,:]
    print N.shape(chisq_array)
    chisq = N.sum(chisq_array)
    n = 1/((2*(N.pi)*(col_sigma_gal**2))**0.5)
    norm = N.sum(N.log(n))
    print 'norm', norm
    sum_chisq = N.sum(chisq_array/2.0, axis=0)
    print 'sum_chisq', sum_chisq
    log_p_d_given_theta = norm  - sum_chisq
    #log_p_d_given_theta = - sum_chisq
    print 'log(P(d|theta))', log_p_d_given_theta
    print 'sum P(d|theta)', N.sum(N.exp(log_p_d_given_theta)) 
    return log_p_d_given_theta

def log_p_d_given_theta(col_gal, col_sigma_gal, col_pred):
    ### col_pred should be a 2D array of dim len(tq) x len(tau) and contain a predicted colour at all of these combinations of parameters for one colour (either nuv-u or u-r)
    ### col_gal and col_sigma_gal should be an array of galaxy colours (either nuv-u or u-r)
    col_gal = col_gal.reshape(len(col_gal), 1)
    col_sigma_gal = col_sigma_gal.reshape(len(col_sigma_gal), 1)
    col_gals = N.split(col_gal, len(col_gal))
    col_sigma_gals = N.array(N.split(col_sigma_gal, len(col_sigma_gal)))
    chisq_array = ((col_gals - col_pred[:,:,0,:])**2)/(col_sigma_gals**2)
    print col_gals - col_pred[:,:,0,:]
    print ((col_gals - col_pred[:,:,0,:])**2)
    print 'chisq', chisq_array[0,:,:]
    print N.shape(chisq_array)
    log_p_d_given_theta = N.sum(chisq_array/2.0, axis=0)
    print 'log(P(d|theta))', log_p_d_given_theta
    print 'sum P(d|theta)', N.sum(N.exp(log_p_d_given_theta))
    return log_p_d_given_theta

def log_p_theta(theta, w):
    mu_tq, mu_tau, sig_tq, sig_tau = w[0], w[1], w[2], w[3]
    theta_mu = theta - N.array([mu_tq,mu_tau]).reshape(2,1)
    C = N.array([1.0/(sig_tq**2), 0.0, 0.0, 1.0/(sig_tau**2)]).reshape(2,2)
    chisq = (theta_mu.T).dot(linalg.inv(C).dot(theta_mu))
    #Z = ((2*N.pi)**(len(C)/2))*((linalg.det(C))**0.5)
    log_p_theta = - (chisq/2.0)
    return log_p_theta

def predict_c_one(theta, age):
    # Time, tq and tau are in units of Gyrs
    print theta
    t = N.linspace(0,14.0,100)
    tq, tau = theta
    sfr = expsfh(tau, tq, t)
    dir ='/Users/becky/Projects/Green-Valley-Project/bc03/models/Padova1994/chabrier/ASCII/'
    model = 'extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
    data = N.loadtxt(dir+model)
    nuv_u = N.zeros_like(sfr)
    u_r = N.zeros_like(sfr)
    total_flux = assign_fluxes.assign_total_flux_numpy(data[0,1:], data[1:,0], data[1:,1:], t*1E9, sfr)
    nuv_u, u_r = get_colours_numpy(t*1E9, total_flux, data)
    nuv_u_age = N.interp(age, t, nuv_u)
    u_r_age = N.interp(age, t, u_r)
    return nuv_u_age, u_r_age

def lnlike_one(theta, c, sigma_c, age):
    tq, tau = theta
    pred_nuv, pred = predict_c_one(theta, age)
    inv_sigma = 1./(sigma_c**2*N.pi*2)**0.5
    return N.log(inv_sigma)-0.5*((c-pred)**2/sigma_c**2)

def lnprior(w, theta):
    mu_tq, mu_tau, sig_tq, sig_tau = w
    tq, tau = theta
    if tq > 0.0 and tq < cosmo.age(0.0) and tau > 0.0 and tau < 3.0:
        ln_tq = -N.log((2*N.pi*sig_tq**2)**0.5) - 0.5*((tq-mu_tq)**2/sig_tq**2)
        ln_tau = -N.log((2*N.pi*sig_tau**2)**0.5) - 0.5*((tau-mu_tau)**2/sig_tau**2)
        return ln_tq + ln_tau
    else:
        return -N.inf

def lnprob(theta, w, c, sigma_c, age):
    lp = lnprior(w, theta)
    if not N.isfinite(lp):
        return -N.inf
    return lp + lnlike_one(theta, c, sigma_c, age)

def sample(ndim, nwalkers, c, sigma_c, age, w):
    pos = [w[0:2] + 1E-3*N.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[w, c, sigma_c, age])
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:,:,:].reshape((-1,ndim))
    fig = triangle.corner(samples, labels=[r'$t$', '$\tau$'])
    fig.savefig('triangle_t_tau.png')
    return samples, fig
    


filters = idlsave.read('/Users/becky/Projects/Green-Valley-Project/Kevin_IDL/ugriz.sav')
fuvwave= filters.ugriz.fuvwave[0]
fuvtrans = filters.ugriz.fuvtrans[0]
nuvwave= filters.ugriz.nuvwave[0]
nuvtrans = filters.ugriz.nuvtrans[0]
uwave= filters.ugriz.uwave[0]
utrans = filters.ugriz.utrans[0]
gwave= filters.ugriz.gwave[0]
gtrans = filters.ugriz.gtrans[0]
rwave= filters.ugriz.rwave[0]
rtrans = filters.ugriz.rtrans[0]
iwave= filters.ugriz.iwave[0]
itrans = filters.ugriz.itrans[0]
zwave= filters.ugriz.zwave[0]
ztrans = filters.ugriz.ztrans[0]
vwave= filters.ugriz.vwave[0]
vtrans = filters.ugriz.vtrans[0]
jwave= filters.ugriz.jwave[0]
jtrans = filters.ugriz.jtrans[0]
hwave= filters.ugriz.hwave[0]
htrans = filters.ugriz.htrans[0]
kwave= filters.ugriz.kwave[0]
ktrans = filters.ugriz.ktrans[0]

def get_colours_numpy(time_steps, sfh, data):
    #fuvmag = get_mag_numpy(time_steps, sfh, fuvwave, fuvtrans, data)
    nuvmag = get_mag_numpy(time_steps, sfh, nuvwave, nuvtrans, data)
    umag = get_mag_numpy(time_steps, sfh, uwave, utrans, data)
    #gmag = get_mag_numpy(time_steps, sfh, gwave, gtrans, data)
    rmag = get_mag_numpy(time_steps, sfh, rwave, rtrans, data)
    #imag = get_mag_numpy(time_steps, sfh, iwave, itrans, data)
    #zmag = get_mag_numpy(time_steps, sfh, zwave, ztrans, data)
    #jmag = get_mag_numpy(time_steps, sfh, jwave, jtrans, data)
    #hmag = get_mag_numpy(time_steps, sfh, hwave, htrans, data)
    #kmag = get_mag_numpy(time_steps, sfh, kwave, ktrans, data)
    nuv_u = nuvmag - umag
    u_r = umag - rmag
    return nuv_u, u_r

def get_mag(total_flux, wave, trans, data):
    mag = N.zeros(len(total_flux))
    for n in range(len(total_flux)):
        mag[n] = assign_fluxes.calculate_AB_mag(data[1:,0], total_flux[n,:], wave, trans)
    return mag

def get_mag_numpy(time_steps, total_flux, wave, trans, data):
    mag = assign_fluxes.calculate_AB_mag_numpy(time_steps, data[1:,0], total_flux, wave, trans)
    return mag
