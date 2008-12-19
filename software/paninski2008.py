#!/usr/bin/env python
#coding=utf-8

import numpy as np
import pylab as plt

def NLIF(v, step, dt, *params):
    g, sigma, I = params  
    return (-g*v + I[step])*dt + sigma*np.sqrt(dt)*np.random.randn()

def FirstPassageMC(func, params, V_thr=1, dt=0.1, max_t=100,
        n_trials=1E4):
    passage_time = []
    time = np.arange(0, max_t, dt)
    for i in xrange(n_trials):
        spike = False
        step, V = 0, 0.
        while (not spike) and (step*dt < max_t):
            dV = func(V, step, dt, *params)
            V = V + dV
            step+=1
            if V>V_thr:
                spike = True
        if spike:
            passage_time.append(step*dt)
    n, bins = np.histogram(passage_time, time)
    return time, n/(1.*n_trials)

def cumsum2d(x):
    """Calculates a matrix of cummulated sums:
       A_{s,t}=\sum_{i=s}^t x_i
    """
    aux = x.cumsum()
    y = np.cumsum(np.tril(np.repeat(aux[:,np.newaxis],len(x), axis=1)),0)
    return y

def cumprod2d(x):
    """Calculates a matrix of cummulated products"""
    n = len(x)
    #aux1 = np.r_[1, x]
    aux = np.repeat(x[:,np.newaxis],n, axis=1)
    aux[np.tri(n, n, -1)==0] = 1 
    y = np.cumprod(aux,0)
    return np.tril(y)

def FirstPassageInt(g, sigma, I, V_thr=1, V_reset=0,dt=0.1, max_t=1000):
    """Calculate first passage time distribution of nLIF using
    integral method 
    (Paninski et al., J Comp Neurosci (2008), 24:69-79)"""

    #1.calculate mean and variance
    u2 = np.exp(-2*g*dt)
    v2 = (1-u2)/(2*g)
    v2[g==0] = dt
    u1 = np.exp(-g*dt)
    v1 = (1-u1)/g
    v1[g==0] = dt

    sigma_sq = np.cumsum(np.vstack((np.zeros((1,len(v2))),v2*cumprod2d(u2))),0)[:-1,:]

    sigma_sq[sigma_sq==0]=1

    mu1 = cumprod2d(u1) 
    mu2 = np.cumsum(np.vstack((np.zeros((1,len(v1))),v1*I*cumprod2d(u1))),0)[:-1,:]
    
    #2. calculate gaussians
    
    def _get_gaussian(y, x):
        mu = mu1*x+mu2
        gaus = 1./np.sqrt(sigma_sq*2*np.pi)*np.exp(-(y-mu)**2/(2*sigma_sq))
        #gaus[isnan(gaus)]=0
        return np.tril(gaus)

    #3. fill in the matrices
    
    A = _get_gaussian(V_thr, V_thr)
    b = _get_gaussian(V_thr, V_reset)[:,0]

    #4. solve the system of linear equations

    p = np.linalg.solve(A,b)

    return p


if __name__ == '__main__':
    #max_t = 5
    #dt = 1
    #g = 1
    #sigma = 1
    #I = 0.5*np.random.randn(int(max_t)/dt)
    #I = 1*np.ones(int(max_t)/dt)
    
    max_t = 50
    dt = 0.5
    g = 0.1
    sigma = 0.1
    #I = 0.5*np.random.randn(int(max_t)/dt)
    I = 0.1*np.ones(np.ceil(max_t/dt))

    p = FirstPassageInt(g*np.ones(len(I)), sigma, I, dt=dt)
    plt.plot(p)

    #t, p_t = FirstPassageMC(NLIF, (g, sigma, I))
    t, p_t = FirstPassageMC(NLIF, (g, sigma, I), V_thr=1,
            dt=dt, max_t=max_t, n_trials=1E4)
    #plt.subplot(211)
    plt.plot(p_t)
    #plt.subplot(212)
    #plt.plot(t, I)
    plt.show()