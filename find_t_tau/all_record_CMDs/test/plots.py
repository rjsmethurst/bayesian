import numpy as N
import pylab as P
import triangle

font = {'family':'serif', 'size':12}
P.rc('font', **font)
P.rc('xtick', labelsize='small')
P.rc('ytick', labelsize='small')

# Load saved text files from bayesian emcee run
times = N.genfromtxt('times.txt', delimiter=',')
print N.shape(times)
urs = N.genfromtxt('opticals.txt', delimiter=',')
print N.shape(urs)
nuvs = N.genfromtxt('nuvs.txt', delimiter=',')
print N.shape(nuvs)
rmags = N.genfromtxt('rmags.txt', delimiter=',')
print N.shape(rmags)
l = N.genfromtxt('likelihoods.txt', delimiter=',')
print N.shape(l)
# find row index of maximum likelihood for each galaxy
idx = N.argmax(l[:,2:], axis=0)

# associate each galaxy with corresponding values for colours, magnitudes etc.
lk = l[:,2:][idx, range(len(idx))] # list of likelihoods len no of galaxies
ts = l[:,0][idx] # list of t_qs len no of galaxies
print N.shape(ts)
taus = l[:,1][idx] # list of tau_qs len no of galaxies
urs = urs[:,2:][idx] # array of urs dim (no of galaxies [y], no of time steps [x])
print N.shape(urs)
nuvs = nuvs[:,2:][idx]  # array of nuvs dim (no of galaxies [y], no of time steps [x])
rmags=rmags[:,2:][idx]  # array of rmags dim (no of galaxies [y], no of time steps [x])

# define past times at which want to investigate CMD
past = N.array([3, 5, 7, 9, 11])

print N.shape(times[2:])

#furs = interpolate.interp2d(times[2:], ts, urs)
#urs_past = furs(past, ts) # array of urs at past times dim (no of galaxies [y], no of past times [x])
#
#fnuvs = interpolate.interp2d(times[2:], ts, nuvs)
#nuvs_past = fuvs(past, ts) # array of nuvs at past times dim (no of galaxies [y], no of past times [x])
#
#frmags = interpolate.interp2d(times[2:], ts, rmags)
#rmags_past = frmags(past, ts) # array of rmags at past times dim (no of galaxies [y], no of past times [x])

urs_past = N.zeros((len(ts), len(past)))
for n in range(len(ts)):
    urs_past[n,:] = N.interp(past, times[2:], urs[n,:])

nuvs_past = N.zeros((len(ts), len(past)))
for n in range(len(ts)):
    nuvs_past[n,:] = N.interp(past, times[2:], nuvs[n,:])

rmags_past = N.zeros((len(ts), len(past)))
for n in range(len(ts)):
    rmags_past[n,:] = N.interp(past, times[2:], rmags[n,:])


fig = P.figure(figsize=(24,5))
for n in range(len(past)):
    ax = P.subplot(1, len(past), n)
    triangle.hist2d(rmags_past[:,n], urs_past[:,n], ax=ax)
    ax.set_xlabel(r'predicted $M_r$')
    #ax.set_ylabel(r'predicted $u-r$ colour')
    ax.set_xlim((-18, -25))
    ax.set_ylim((-3, 4))
ax = P.subplot(1, 5, 1)
ax.set_ylabel(r'predicted $u-r$ colour')
P.tight_layout()
fig.savefig('cmd_evo.pdf')

fig = P.figure(figsize=(24,5))
for n in range(len(past)):
    ax = P.subplot(1, len(past), n)
    triangle.hist2d(urs_past[:,n], nuvs_past[:,n], ax=ax)
    #ax.set_ylabel(r'predicted $NUV-u$ colour')
    ax.set_xlabel(r'predicted $u-r$ colour')
    ax.set_xlim((-3,4))
    ax.set_ylim((-3,4))
ax = P.subplot(1, 5, 1)
ax.set_ylabel(r'predicted $NUV-u$ colour')
P.tight_layout()
fig.savefig('c_c_evo.pdf')
    
    
    




