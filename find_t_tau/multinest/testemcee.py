import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.optimize as opt
import time

plt.ion()

'''Fit a line using the MCMC Hammer - just a quick demonstration'''

def model(m,b,x):
	'''In this case, a straight line'''
	return x*m+b

def lnprob(walker,x,y,error):
	chi2 = np.sum((model(walker[0],walker[1],x)-y)**2/error**2)
	return -chi2/2.

xx = np.linspace(0,30,100)
mm = 5.0
bb = -3.0
noise = 10.0

yy = mm*xx+bb + np.random.randn(len(xx))*noise

#initialise
ndim, nwalkers = 2, 100
nsteps = 1000
ivar = 1. / np.random.rand(ndim)
p0 = [np.random.rand(ndim) for i in range(nwalkers)]
tic = time.time()
#run!
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[xx,yy,noise])

# burn in
pos,prob,state = sampler.run_mcmc(p0, 100)
sampler.reset()

# restart
sampler.run_mcmc(pos,nsteps)
toc = time.time()
print 'Time elapsed:',toc-tic,'s'

#get histograms
ms = sampler.flatchain[:,0]
bs = sampler.flatchain[:,1]

print 'Gradient =',np.mean(ms),'pm',np.std(ms)
print 'Intercept =',np.mean(bs),'pm',np.std(bs)


'''--------------------------
Generate a video of the chain
--------------------------'''

for k in range(0,nsteps/5):
 	print 'generating plot'
 	plt.clf()
 	start = k*5*nwalkers
 	plt.scatter(ms[start:start+nwalkers],bs[start:start+nwalkers])
 	plt.title('MCMC Hammer Fit to a Simple Function')
 	plt.xlabel('Gradient')
 	plt.ylabel('y-intercept')
 	plt.axis([ms.min(),ms.max(),bs.min(),bs.max()])
 	plt.draw()
 	plt.show()

plt.clf()
plt.scatter(ms,bs)
plt.title('MCMC Hammer Fit to a Simple Function')
plt.xlabel('Gradient')
plt.ylabel('y-intercept')
plt.axis([ms.min(),ms.max(),bs.min(),bs.max()])
plt.draw()
plt.show()

plt.figure(2)
plt.clf()
plt.errorbar(xx,yy,yerr=noise,fmt=' ')
plt.plot(xx,model(np.mean(ms),np.mean(bs),xx))
plt.legend(['Data','Best Fit'])
plt.title('MCMC Hammer Fit to a Simple Function')
plt.xlabel('Ordinate')
plt.ylabel('Abscissa')
plt.show()