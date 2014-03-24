import numpy as N
import pylab as P
import scipy as S
import pyfits as F
from t_tau_func import *
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0 = 71.0, Om0 = 0.26)
av_z = 0.076401
age = cosmo.age(av_z)


font = {'family':'serif', 'size':12}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')

dir ='/Users/becky/Projects/Green-Valley-Project/bc03/models/Padova1994/chabrier/ASCII/'
model = 'extracted_bc2003_lr_m62_chab_ssp.ised_ASCII'
data = N.loadtxt(dir+model)

tq = N.linspace(0.0, 13.6, 10)
tau = N.linspace(3.0, 0.001, 10)

ur = N.zeros((len(tq),len(tau)))
nuv = N.zeros_like(ur)
for n in range(len(tq)):
    for m in range(len(tau)):
        nuv[n,m], ur[n,m] = predict_c_one([tq[n], tau[m]], age)

N.save('ur.npy', ur)
N.save('nuvu.npy', nuv)

print 'u-r',ur
print 'nuv-u', nuv

P.figure()
#P.imshow(ur, origin ='lower', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)], aspect='auto', interpolation='nearest')
P.imshow(ur, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$')
cbar = P.colorbar()
cbar.set_label(r'predicted $u-r$ colour')
P.show()

P.figure()
P.imshow(nuv, origin='lower', aspect='auto', extent=[N.min(tq), N.max(tq), N.min(tau), N.max(tau)])
P.xlabel(r'$t_{quench} (Gyr)$')
P.ylabel(r'$\tau$')
cbar= P.colorbar()
cbar.set_label(r'predicted $NUV-u$ colour')
P.show()
        