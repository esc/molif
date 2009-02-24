#!/usr/bin/env python
#coding=utf-8

# Copyright 2008 Valentin 'esc' Haenel
#
# This program is free software. It comes without any warranty, to
# the extent permitted by applicable law. You can redistribute it
# and/or modify it under the terms of the Do What The Fuck You Want
# To Public License, Version 2, as published by Sam Hocevar. See
# http://sam.zoy.org/wtfpl/COPYING for more details. 

from pylab import imshow, figure, plot, show, colorbar, hold, \
xlabel, ylabel, legend, subplot, title
from numpy import flipud, cumsum
from matplotlib import cm
from scipy import interpolate

from density import *
from model import *
from montecarlo import *

def plot_pde_P_vt_fpt():
    """ calculate and plot density evolution and and first passage time. """
    P_vt, fpt = compute_pde_fpt()
    imshow(flipud(P_vt))
    colorbar()
    figure()
    plot(fpt)
    show()

def plot_mc_P_vt_fpt():
    """ compute and plot monte carlo density evolution and first
    passage time"""

    P_vt, traces, fpt =  mc_P_vt_fpt()

    imshow(P_vt)
    colorbar()

    figure()
    hold(True)
    for t in traces:
        plot(t)
    xlabel('time')
    ylabel('V(t)')

    figure()
    plot(fpt)
    show()


def plot_mc_fpt():
    """ compute and plot the monte carlo first passage time """
    fpt = mc_fpt(reps=5000)
    plot(fpt)
    show()

def compare_pde_mc_P_vt():


    pde_P_vt_val, pde_fpt_val = compute_pde_fpt()
    mc_P_vt_val, traces, mc_fpt_val =  mc_P_vt_fpt()

    subplot(2,1,1)
    imshow(flipud(pde_P_vt_val[375:,:250]),extent=(0.0,250.0,0.0,1.0), aspect='auto')
    colorbar()
    xlabel('time')
    ylabel('P(V,t)')
    title('PDE')

    subplot(2,1,2)
    imshow(mc_P_vt_val[:125,:250],extent=(0.0,250.0,0.0,1.0), aspect='auto')
    colorbar()
    xlabel('time')
    ylabel('P(V,t)')
    title('Monte Carlo')

    show()
 
def compare_pde_mc_fpt():
    """ compare the partial differental equation and monte carlo first
    passage time """
    lif = lif_setup()
    n_repetitions = 5000

    mc_time, mc_fpt, spike_times = compute_mc_fpt(reps=n_repetitions)
    pde_time, pde_fpt = compute_pde_fpt()
    #integral_time, integral_fpt = compute_integral_fpt()
    
    mc_fpt = mc_fpt[:len(pde_fpt)]
    mc_time = mc_time[:len(pde_fpt)]

    plot(mc_time, mc_fpt,'rx',label='Monte Carlo')
    plot(pde_time, pde_fpt,'gx',label='PDE Solver')
    #plot(integral_time, integral_fpt/5.,'b',label='Integral')
    xlabel('time')
    ylabel('fpt')
    legend()

    cdf = interpolate.interp1d(pde_time, cumsum(pde_fpt)*lif.dt)
    D, p = stats.kstest(spike_times, cdf)
    print "K-S test D value: " , D
    print "K-S test p value: " , p


    [n, bins] = histogram(spike_times, 100)
    mc_cdf = cumsum(n)*1./n_repetitions
    figure()
    
    plot(bins[:-1], mc_cdf,'b')
    plot(pde_time, cdf(pde_time),'g')
    show()



if __name__ == '__main__':
    #compare_pde_mc_fpt()
    plot_mc_P_vt_fpt()
