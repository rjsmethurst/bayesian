import numpy as N
import scipy as S
import pylab as P
import pyfits as F
import idlsave
from scipy.integrate import simps
import scipy.spatial.distance as SD
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
import time

cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)

########################################################################################
# The data files .ised_ASCII contain the extracted bc03 models and have a 0 in the origin at [0,0]. The first row contains
# the model ages (from the second column) - data[0,1:]. The first column contains the model lambda values (from the second
# row) - data[1:,0]. The remaining data[1:,1:] are the flux values at each of the ages (columns, x) and lambda (rows, y)
# values.
########################################################################################


dir ='/Users/becky/Projects/Green-Valley-Project/bc03/models/Padova1994/chabrier/ASCII/'
model = 'extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
data = N.loadtxt(dir+model)

# Function which given a tau and a tq calculates the sfr at all times
def expsfh(tau, tq, time):
    a = time.searchsorted(tq)
    sfr = N.ones(len(time))*100
    sfr[a:] = 100*N.exp(-(time[a:]-tq)/tau)
    return sfr

# predict the colour of a galaxy of a given age given a sf model of tau and tq
def predict_c_one(theta, age):
    # Time, tq and tau are in units of Gyrs
    time = N.arange(0, 0.01, 0.003)
    t = N.linspace(0,14.0,200)
    t = N.append(time, t[1:])
    tq, tau = theta
    sfr = expsfh(tau, tq, t)
    nuv_u = N.zeros_like(sfr)
    u_r = N.zeros_like(sfr)
    # Work out total flux at each time given the sfh model of tau and tq (calls assign_fluxes function)
    total_flux = assign_fluxes.assign_total_flux(data[0,1:], data[1:,0], data[1:,1:], t*1E9, sfr)
    # Calculate fluxes from the flux at all times then interpolate to get one colour for the age you are observing the galaxy at - if many galaxies are being observed, this also works with an array of ages to give back an array of colours
    nuv_u, u_r = get_colours(t*1E9, total_flux, data)
    nuv_u_age = N.interp(age, t, nuv_u)
    u_r_age = N.interp(age, t, u_r)
    return nuv_u_age, u_r_age

#Calculate colours and magnitudes for functions above
def get_colours(time_steps, sfh, data):
    nuvmag = get_mag(time_steps, sfh, nuvwave, nuvtrans, data)
    umag = get_mag(time_steps, sfh, uwave, utrans, data)
    rmag = get_mag(time_steps, sfh, rwave, rtrans, data)
    nuv_u = nuvmag - umag
    u_r = umag - rmag
    return nuv_u, u_r

def get_mag(time_steps, total_flux, wave, trans, data):
    mag = assign_fluxes.calculate_AB_mag(time_steps, data[1:,0],total_flux, wave, trans)
    return mag

# Function for likelihood of model given all galaxies
def lnlike_one(theta, ur, sigma_ur, nuvu, sigma_nuvu, age):
    tq, tau = theta
    pred_nuvu, pred_ur = predict_c_one(theta, age)
    #inv_sigma_ur = 1./((sigma_ur**2)*(2*N.pi))**0.5
    #inv_sigma_nuvu = 1./((sigma_nuvu**2)*(2*N.pi))**0.5
    #return N.log(inv_sigma)-0.5*((ur-pred_ur)**2/sigma_ur**2) #- 0.5*((nuvu-pred_nuvu)**2/sigma_nuvu**2)
    return -0.5*((ur-pred_ur)**2/sigma_ur**2)-0.5*((nuvu-pred_nuvu)**2/sigma_nuvu**2)

# Function which includes GZ likelihoods and sums across all galaxies to return one value for a given set of theta 
def lnlike(theta, ur, sigma_ur, nuvu, sigma_nuvu, age, pd, ps):
    ts, taus, td, taud = theta
    d = lnlike_one([td, taud], ur, sigma_ur, nuvu, sigma_nuvu, age)
    s = lnlike_one([ts, taus], ur, sigma_ur, nuvu, sigma_nuvu, age)
    D = N.log(pd) + d
    S = N.log(ps) + s
    return N.sum(N.logaddexp(D, S))

# Prior likelihood on theta values given the inital w values assumed for the mean and stdev
def lnprior(w, theta):
    mu_tqs, mu_taus, mu_tqd, mu_taud, sig_tqs, sig_taus, sig_tqd, sig_taud = w
    ts, taus, td, taud = theta
    if 0.0 < ts < 13.807108309208775 and 0.0 < taus < 3.0 and 0.0 < td < 13.807108309208775 and 0.0 < taud < 3.0:
        ln_tqs = - 0.5*((ts-mu_tqs)**2/sig_tqs**2) #- N.log((2*N.pi*sig_tqs**2)**0.5)
        ln_taus = - 0.5*((taus-mu_taus)**2/sig_taus**2) #-N.log((2*N.pi*sig_taus**2)**0.5)
        ln_tqd = - 0.5*((td-mu_tqd)**2/sig_tqd**2) #-N.log((2*N.pi*sig_tqd**2)**0.5) 
        ln_taud = - 0.5*((taud-mu_taud)**2/sig_taud**2) #-N.log((2*N.pi*sig_taud**2)**0.5)
        #print 'prior', ln_tqs + ln_taus + ln_tqd + ln_taud
        return ln_tqs + ln_taus + ln_tqd + ln_taud
    else:
        return -N.inf

# Overall likelihood function combining prior and model
def lnprob(theta, w, ur, sigma_ur, nuvu, sigma_nuvu, age, pd, ps):
    lp = lnprior(w, theta)
    if not N.isfinite(lp):
        return -N.inf
    return lp + lnlike(theta, ur, sigma_ur, nuvu, sigma_nuvu, age, pd, ps)


def sample(ndim, nwalkers, w, ur, sigma_ur, nuvu, sigma_nuvu, age, pd, ps):
    if len(age) != len(ur):
        raise SystemExit('Number of ages does not coincide with number of galaxies...')
    pos = [[8.0, 1.25, 8.0, 1.25] + 1e-4*N.random.randn(ndim) for i in range(nwalkers)]
    n = 0
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(w, ur, sigma_ur, nuvu, sigma_nuvu, age, pd, ps))
    sampler.run_mcmc(pos, 500)
    samples = sampler.chain[:,50:,:].reshape((-1,ndim))
    samples_save = '/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/blackbox/samples_all_'+str(len(samples))+'_'+str(len(age))+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.npy'
    N.save(samples_save, samples)
    fig = triangle.corner(samples, labels=[r'$ t_{smooth} $', r'$ \tau_{smooth} $', r'$ t_{disc} $', r'$ \tau_{disc}$'])
    fig.savefig('triangle_t_tau_all_'+str(len(samples))+'_'+str(len(age))+'_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf')
    return samples, fig


#Define function to plot the walker positions as a function of the step
def walker_plot(samples, nwalkers, limit):
    s = samples.reshape(nwalkers, -1, 4)
    s = s[:,:limit, :]
    fig = P.figure(figsize=(8,10))
    ax1 = P.subplot(4,1,1)
    ax2 = P.subplot(4,1,2)
    ax3 = P.subplot(4,1,3)
    ax4 = P.subplot(4,1,4)
    for n in range(len(s)):
        ax1.plot(s[n,:,0], 'k')
        ax2.plot(s[n,:,1], 'k')
        ax3.plot(s[n,:,2], 'k')
        ax4.plot(s[n,:,3], 'k')
    ax1.tick_params(axis='x', labelbottom='off')
    ax2.tick_params(axis='x', labelbottom='off')
    ax3.tick_params(axis='x', labelbottom='off')
    ax4.set_xlabel(r'step number')
    ax1.set_ylabel(r'$t_{smooth}$')
    ax2.set_ylabel(r'$\tau_{smooth}$')
    ax3.set_ylabel(r'$t_{disc}$')
    ax4.set_ylabel(r'$\tau_{disc}$')
    P.subplots_adjust(hspace=0.1)
    save_fig = '/Users/becky/Projects/Green-Valley-Project/bayesian/find_t_tau/blackbox/walkers_steps_all_'+str(time.strftime('%H_%M_%d_%m_%y'))+'.pdf'
    fig.savefig(save_fig)
    return fig


#Load the filters in order to calculate fluxes in each bandpass
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
