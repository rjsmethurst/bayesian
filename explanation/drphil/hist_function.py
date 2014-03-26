import numpy as N
import pylab as P
import pyfits as F

def bin(xmin, xmax, ymin, ymax, xnumbins, ynumbins, colours):
    xbinsep = (xmax-xmin)/xnumbins
    ybinsep = (ymax-ymin)/ynumbins
    xbins = N.arange(xmin, xmax+xbinsep, xbinsep)
    ybins = N.arange(ymin, ymax+ybinsep, ybinsep)
    grid = N.zeros(len(xbins)*len(ybins)).reshape(len(ybins), len(xbins))
    for n in range (len(xbins)-1):
        for m in range(len(ybins)-1):
            grid[m, n] = ((xbins[n] < colours[:,0]) & (colours[:,0] < xbins[n+1]) & (ybins[m] < colours[:,1]) & (colours[:,1] < ybins[m+1])).sum()
    return xmin, xmax, ymin, ymax, xbins, ybins, grid

def xy_bin(xmin, xmax, ymin, ymax, xnumbins, ynumbins, x, y):
    xbinsep = (xmax-xmin)/xnumbins
    ybinsep = (ymax-ymin)/ynumbins
    xbins = N.arange(xmin, xmax+xbinsep, xbinsep)
    ybins = N.arange(ymin, ymax+ybinsep, ybinsep)
    grid = N.zeros(len(xbins)*len(ybins)).reshape(len(ybins), len(xbins))
    for i in range (len(xbins)-1):
        for j in range(len(ybins)-1):
            grid[i, j] = ((xbins[i] < x) & (x < xbins[i+1]) & (ybins[j] < y) & (y < ybins[j+1])).sum()
    return xmin, xmax, ymin, ymax, xbins, ybins, grid

def avg_bin(x, xres, y, yres, z):
    ybins = N.linspace(0, 13.8, yres)
    xbins = N.linspace(0, 3.0, xres)
    X, Y = N.meshgrid(xbins, ybins)
    grid = N.zeros((len(xbins)-1)*(len(ybins)-1)).reshape((len(ybins)-1), (len(xbins)-1))
    for i in range(len(xbins)-1):
        for j in range(len(ybins)-1):
            zlist=[]
            for n in range (len(z)-1):
                if x[n] >= xbins[i] and x[n] <= xbins[i+1] and y[n] >= ybins[j] and y[n] <= ybins[j+1]:
                    zlist = N.append(zlist, z[n])
            grid[len(ybins)-2-j,i] = N.mean(zlist)
    return X, Y, grid

def sd_bin(x, xres, y, yres, z):
    ybins = N.linspace(-2, N.max(y)+0.5, yres)
    xbins = N.linspace(N.min(x)-0.5, 4.5, xres)
    X, Y = N.meshgrid(xbins, ybins)
    grid = N.zeros((len(xbins)-1)*(len(ybins)-1)).reshape((len(ybins)-1), (len(xbins)-1))
    for i in range(len(xbins)-1):
        for j in range(len(ybins)-1):
            zlist=[]
            for n in range (len(z)-1):
                if x[n] >= xbins[i] and x[n] <= xbins[i+1] and y[n] >= ybins[j] and y[n] <= ybins[j+1]:
                    zlist = N.append(zlist, z[n])
            grid[len(ybins)-2-j,i] = N.std(zlist)
    return X, Y, grid


def griddata(x, xres, y, yres, z):
    yi = N.linspace(-2, N.max(y)+0.5, yres)
    xi = N.linspace(N.min(x)-0.5, 4.5, xres)
    X, Y = N.meshgrid(xi,yi)
    Z = P.griddata(x,y,z,xi, yi)
    return X, Y, Z

def probgrid(px,py,x,y):
    extent = (min(x), max(x), min(y), max(y))
    prob = N.outer(py,px)
    return prob, extent

def normaldist(x, u, s):
    if s > 0:
        p_x = N.zeros_like(x)
        for n in range(len(x)):
            p_x[n] = (1/((2*N.pi*(s**2))**0.5))*N.exp(-(((x[n]-u)**2)/(2*(s**2))))
    else:
        p_x=[1.0]
    return p_x

def chisq(di, dk, si):
    chi = N.zeros(len(dk))
    for n in range(len(dk)):
        chi[n] = (((di - dk[n])**2)/(si)**2)
    chisq = chi.sum()
    p_d_theta = (1.0/(2*N.pi*(si**2)))*N.exp(-chisq/2.0)
    log_p_d_theta = p_d_theta
    if log_p_d_theta == N.inf:
        return 1.0
    if log_p_d_theta == -N.inf:
        return 0.0
    else:
        return log_p_d_theta


