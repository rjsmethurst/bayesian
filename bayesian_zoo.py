import scipy as S
import numpy as N
import pylab as P
from matplotlib import pyplot as PY
import idlsave
from hist_function import *
import pyfits as F
from scipy.stats import norm, gumbel_r, gumbel_l
from scipy import linalg

font = {'family':'serif', 'size':12}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')

file = '/Users/becky/Projects/Green-Valley-Project/data/GZ2_all_GALEX_match.fits'
dat = F.open(file)
gz2data = dat[1].data
col = N.zeros(6*len(gz2data)).reshape(len(gz2data), 6)
col[:,0] = gz2data.field('petroMag_u') - gz2data.field('petroMag_r')
col[:,1] = gz2data.field('NUV') - gz2data.field('petroMag_u')
col[:,3] = (gz2data.field('e_NUV') + gz2data.field('petroMagErr_u'))
col[:,2] = (gz2data.field('petroMagErr_r') + gz2data.field('petroMagErr_u'))
col[:,4] = gz2data.field('t01_smooth_or_features_a01_smooth_flag')
col[:,5] = gz2data.field('t01_smooth_or_features_a02_features_or_disk_flag')

non_nan = N.logical_not(N.isnan(col[:,1])).astype(int)
colours = N.compress(non_nan, col, axis=0)

disc = colours[colours[:,5]==1]
smooth = colours[colours[:,4]==1]
inter = colours[colours[:,4]!=1]
intermediate = inter[inter[:,5]!=1]

sfhs = ['exp_25gyr','exp_1gyr', 'exp_250myr', 'exp_1myr', 'truncated', 'constant', 'burst']
sfh_label = ['2.5 Gyr','1 Gyr', '250 Myr', '1 Myr', 'Truncated', 'Constant', 'Burst with Exponential 2.5 Gyr']
c = ['g', 'r', 'm', 'b', 'k', 'y', 'c']
metal =[62]

for m in range(len(sfhs)-3):
    for n in range(len(metal)):
        globals()['model_sfh_'+sfhs[m]+'_'+str(metal[n])] = N.loadtxt('time_colour_array_'+str(metal[n])+'_'+sfhs[m])

#P.figure()
#for n in range(len(sfhs)-3):
#    P.plot(globals()['model_sfh_'+sfhs[n]+'_62'][:,2], globals()['model_sfh_'+sfhs[n]+'_62'][:,1], c[n]+'o-', markevery=10)
#    P.plot(globals()['model_sfh_'+sfhs[n]+'_62'][model_sfh_exp_25gyr_62[:,0]==9.0][:,2], globals()['model_sfh_'+sfhs[n]+'_62'][model_sfh_exp_25gyr_62[:,0]==9.0][:,1], c[n]+'D', markersize=10)
#P.show()

def logptheta(theta, w):
    mu_tq, mu_tau, sig_tq, sig_tau = w[0], w[1], w[2], w[3]
    theta_mu = theta - N.array([mu_tq,mu_tau]).reshape(2,1)
    C = N.array([1.0/sig_tq, 0.0, 0.0, 1.0/sig_tau]).reshape(2,2)
    print linalg.inv(C)
    print linalg.inv(C).dot(theta_mu)
    chisq = (theta_mu.T).dot(linalg.inv(C).dot(theta_mu))
    print chisq
    Z = ((2*N.pi)**(len(C)/2))*((linalg.det(C))**0.5)
    print Z
    log_p_theta = -N.log(Z) - (chisq/2.0)
    return log_p_theta

tau = N.arange(0.0, 2.5, 0.05)
tq = N.arange(5.0, 13.0, 0.1)
ext = (N.min(tau), N.max(tau), N.min(tq), N.max(tq))
print ext


w = [9.0, 1.0, 0.5, 0.1]
prob = N.zeros(len(tau)*len(tq)).reshape(len(tq), len(tau))
for n in range(len(tau)):
    for m in range(len(tq)):
        theta = N.array([tq[m], tau[n]]).reshape(2,1)
        prob[m,n] = logptheta(theta, w)


#TAU, TQ = N.meshgrid(tau, tq)
#
dist_tq = norm(loc=w[0], scale=w[2])
dist_tau = norm(loc=w[1], scale = w[3])
#
pdf_tau = dist_tau.pdf(tau)
pdf_tq = dist_tq.pdf(tq)
#
#prob, ext = probgrid(pdf_tau, pdf_tq, tau, tq)

lev = [0.1, 0.3, 0.5, 0.7, 0.9]
levels = lev*N.max(prob)
#
#
ax1 = PY.subplot(2,2,2)
cont = ax1.contour(prob, extent=ext, cmap=PY.cm.binary)
ax1.set_xlabel(r'$\tau  (Gyr)$')
ax1.set_ylabel(r'$t_{quench}  (Gyr)$')
ax2 = PY.subplot(2,2,1)
ptq = P.plot(pdf_tq, tq, color='k')
ptq = PY.gca()
ptq.invert_xaxis()
ax2.tick_params(axis='y', labelright='off', labelleft='off')
ax2.set_xlabel('Probability Density')
ax3 = PY.subplot(2,2,4)
ptau = PY.plot(tau, pdf_tau, color='k')
ptau = PY.gca()
ptau.invert_yaxis()
ax3.tick_params(axis='x', labeltop='off', labelbottom='off')
ax3.set_ylabel('Probability Density')
PY.subplots_adjust(wspace = 0.25, hspace =0.25)
PY.show()


p_d_theta = N.zeros(10000)
si = N.linspace(0.5, 0.5, len(p_d_theta))
for n in range(len(p_d_theta)):
    p_d_theta[n] = chisq(disc[n,0], model_sfh_exp_25gyr_62[model_sfh_exp_25gyr_62[:,0]==9.0][:,2], disc[n,2])

#P.figure()
#P.plot(colours[0:len(p_d_theta),0], p_d_theta, 'kx')
#P.xlabel(r'$u-r$')
#P.ylabel(r'$P(d_{k}|\theta_{k}, w, H)$')
#P.ylim(ymax=1.0)
#P.xlim(xmin=0.0, xmax=6.0)
#P.show()
