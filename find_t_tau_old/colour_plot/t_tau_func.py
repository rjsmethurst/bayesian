import numpy as N
import scipy as S
import pylab as P
import pyfits as F
import idlsave
from astropy.cosmology import FlatLambdaCDM
import assign_fluxes
import pyfits as F
import time

cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)

#########################################################################################################################
# The data files .ised_ASCII contain the extracted bc03 models and have a 0 in the origin at [0,0]. The first row contains
# the model ages (from the second column) - data[0,1:]. The first column contains the model lambda values (from the second
# row) - data[1:,0]. The remaining data[1:,1:] are the flux values at each of the ages (columns, x) and lambda (rows, y)
# values.
#########################################################################################################################

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
    t = N.linspace(0,14.0,100)
    tq, tau = theta
    sfr = expsfh(tau, tq, t)
    nuv_u = N.zeros_like(sfr)
    u_r = N.zeros_like(sfr)
    # Work out total flux at each time given the sfh model of tau and tq (calls assign_fluxes function)
    # The total_flux array has time_steps as the rows (y) and model_lambda along the columns(x)
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
