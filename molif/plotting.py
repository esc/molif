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

from density import *
from model import *
from montecarlo import *

def pde_density_and_fpt():
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
    for t in traces:  plot(t)

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
    imshow(flipud(pde_P_vt_val[375:,:250]), cmap=cm.Greys_r)
    xlabel('time')
    ylabel('P(V,t)')
    title('PDE')

    subplot(2,1,2)
    imshow(mc_P_vt_val[70:,:250], cmap=cm.Greys_r)
    xlabel('time')
    ylabel('P(V,t)')
    title('Monte Carlo')

    show()
 
def compare_pde_mc_fpt():
    """ compare the partial differental equation and monte carlo first
    passage time """

    mc_fpt_val = mc_fpt(reps=5000)
    P_vt, pde_fpt = compute_pde_fpt()

    D,p = stats.ks_2samp(mc_fpt_val,pde_fpt)
    print "K-S test D value: " , D
    print "K-S test p value: " , p
    plot(mc_fpt_val,'r',label='Monte Carlo')
    plot(pde_fpt,'g',label='PDE Solver')
    xlabel('time')
    ylabel('fpt')
    legend()


    figure()

    plot(cumsum(mc_fpt_val),'r')
    plot(cumsum(pde_fpt),'g')
    show()



if __name__ == '__main__':
    compare_pde_mc_P_vt()
