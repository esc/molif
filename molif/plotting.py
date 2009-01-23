#!/usr/bin/env python
#coding=utf-8

# Copyright 2008 Valentin 'esc' Haenel
#
# This program is free software. It comes without any warranty, to
# the extent permitted by applicable law. You can redistribute it
# and/or modify it under the terms of the Do What The Fuck You Want
# To Public License, Version 2, as published by Sam Hocevar. See
# http://sam.zoy.org/wtfpl/COPYING for more details. 

from pylab import imshow, figure, plot, show, colorbar
from numpy import flipud, cumsum

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

def plot_monte_carlo_fpt():

    mc_fpt = monte_carlo_fpt(reps=5000)
    plot(mc_fpt)
    show()


def compare_pde_mc_fpt():
    """ compare the partial differental equation and monte carlo first
    passage time """

    mc_fpt = monte_carlo_fpt(reps=5000)
    P_vt, pde_fpt = compute_pde_fpt()

    D,p = stats.ks_2samp(mc_fpt,pde_fpt)
    print "K-S test D value: " , D
    print "K-S test p value: " , p
    plot(mc_fpt,'r')
    plot(pde_fpt,'g')
    figure()

    plot(cumsum(mc_fpt),'r')
    plot(cumsum(pde_fpt),'g')
    show()




