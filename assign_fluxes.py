import numpy as N
import scipy as S
from scipy.integrate import simps
from scipy import interpolate
import pylab as P

def assign_total_flux(model_ages, model_lambda, model_fluxes, time_steps, sim_SFR):
    #print model_ages
    mask = model_ages[model_ages<4E6]
    model_fluxes[:,0:len(mask)] = 0.0    
    interp_fluxes_sim = N.zeros(len(model_lambda)*len(time_steps)).reshape(len(model_lambda), len(time_steps))
    for n in range(len(model_fluxes)):
        interp_fluxes_sim[n,:] = N.interp(time_steps, model_ages, model_fluxes[n,:])
    fraction_array = N.zeros(len(time_steps)*len(time_steps)).reshape(len(time_steps), len(time_steps))
    for j in range(len(time_steps)): # Over all columns
        fraction_array[0,j] = sim_SFR[j]/sim_SFR[0]
        for i in range(1,len(time_steps)+1): # Over all rows
            if i <= j: # Only for parts of the array where model ages aren't greater than the time_step
                fraction_array[i,j] = fraction_array[i-1, j-1]
            else:
                pass
    mass_array = N.zeros(len(time_steps)*len(time_steps)).reshape(len(time_steps), len(time_steps))
    for j in range(len(time_steps)): # Over all columns
        mass_array[0,0] = sim_SFR[0]
        #mass_array[0,j] = S.integrate.quad(sim_SFR, time_step[j-1:j], dt)
        mass_array[0,j] = sim_SFR[j]*(time_steps[j]-time_steps[j-1])
        for i in range(1,len(time_steps)+1): # Over all rows
            if i <= j: # Only for parts of the array where model ages aren't greater than the time_step
                mass_array[i,j] = mass_array[i-1, j-1]
            else:
                pass
    frac_flux_array = fraction_array*mass_array
    total_flux=N.zeros(len(model_lambda)*len(time_steps)).reshape(len(time_steps), len(model_lambda))
    for n in range(len(model_lambda)):
        flux_array = N.zeros_like(frac_flux_array)
        for i in range(len(time_steps)):
            flux_array[i,:] = frac_flux_array[i,:]*interp_fluxes_sim[n,i]            
        total_flux[:,n] = N.sum(flux_array, axis=0)
    return total_flux

def assign_total_flux_numpy(model_ages, model_lambda, model_fluxes, time_steps, sim_SFR):
        #First mask the ages of the very young stars hidden in birth clouds
    mask = model_ages[model_ages<4E6]
    model_fluxes[:,0:len(mask)] = 0.0
        # Calculate the fluxes at the ages specified by the time steps rather than in the models using numpy/scipy array manipulations rather than a for loop
    f = interpolate.interp2d(model_ages, model_lambda, model_fluxes)
    interp_fluxes_sim = f(time_steps, model_lambda)
        # Produce the array to keep track of the ages of the fractional SFR at each time step
    frac_sfr = sim_SFR/sim_SFR[0]
    fraction_array = S.linalg.toeplitz(frac_sfr, N.zeros_like(frac_sfr)).T
        # Produce the array to keep track of the ages of the mass fraction of stars formed at each time step
    m_array = (sim_SFR.T)*(N.append(1, N.diff(time_steps)))
    mass_array = S.linalg.toeplitz(m_array, N.zeros_like(frac_sfr)).T
        # Produce the array to keep track of the fraction of flux produced at each timestep 
    frac_flux_array = fraction_array*mass_array
        # Calculate the total flux contributed by all of the fractions at each time step by summing across all wavelength values
    flux_lambda = frac_flux_array*(N.split(interp_fluxes_sim.T, len(model_lambda), axis=1))
    total_flux = (N.sum(flux_lambda, axis=1)).T # Array of dimensions (len(timesteps), len(model_lambda))
    return total_flux


def dust_calzetti(sed_lambda, sed_flux, ebmv):
    k = N.zeros_like(sed_lambda)
    for i in range(len(sed_lambda)-1):
        l = sed_lambda/1E4
        if l[i] > 0.63:
            k[i] = 2.659*(-1.857 + 1.040/l[i]) + 4.05
        if l[i] < 0.63:
            k[i] = 2.659*(-2.156 + (1.509/l[i]) - (0.198/(l[i]**2)) + (0.011/(l[i]**3))) + 4.05
        if l[i] > 2.2:
            k[i] = 0
    dust_corr_flux = sed_flux*(10**(-0.4*ebmv*k))
    return dust_corr_flux


def calculate_AB_mag(model_lambda, sim_flux, wave, trans):

    lambda_filter1 = [i for i in model_lambda if i > wave[0] and i < wave[len(wave)-1]]
    lambda_filter2 = N.append(wave[0], lambda_filter1)
    lambda_filter = N.append(lambda_filter2, wave[len(wave)-1])
    
    flux_filter = N.interp(lambda_filter, model_lambda, sim_flux)
    trans_filter = N.interp(lambda_filter, wave, trans)

    top = N.trapz(lambda_filter*flux_filter*trans_filter, lambda_filter)
    bottom = N.trapz(trans_filter/lambda_filter, lambda_filter)

    m_ab = -2.41 - 2.5*N.log10(top/bottom)
    
    return m_ab

def calculate_AB_mag_numpy(time_steps, model_lambda, sim_flux, wave, trans):
    
    lambda_filter1 = [i for i in model_lambda if i > wave[0] and i < wave[len(wave)-1]]
    lambda_filter2 = N.append(wave[0], lambda_filter1)
    lambda_filter = N.append(lambda_filter2, wave[len(wave)-1])

    f = interpolate.interp2d(model_lambda, time_steps, sim_flux)
    flux_filter = (f(lambda_filter, time_steps))
    trans_filter = N.interp(lambda_filter, wave, trans)
    
    top = N.trapz((lambda_filter*flux_filter*trans_filter), lambda_filter, axis=1)
    bottom = N.trapz(trans_filter/lambda_filter, lambda_filter)
    
    m_ab = -2.41 - 2.5*N.log10(top/bottom)
    
    return m_ab
            


