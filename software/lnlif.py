#!/usr/bin/env python
#coding=utf-8

# Copyright 2008 Valentin 'esc' Haenel
#
# This program is free software. It comes without any warranty, to
# the extent permitted by applicable law. You can redistribute it
# and/or modify it under the terms of the Do What The Fuck You Want
# To Public License, Version 2, as published by Sam Hocevar. See
# http://sam.zoy.org/wtfpl/COPYING for more details. 

from numpy import arange, zeros, sqrt, random, diff, cumsum
from scipy import stats, sparse, linsolve, optimize



def plot_three_h():
    """ shows the possibilities that we have with varying h """

    lif = lnlif() # init model
    lif.set_const_input(0.01); # set constant input
    lif.i_stim = lif.stim # setup stimulus
    # lif.set_convolved_input();
    #lif.noise = True
    lif.set_const_h();

    time, potential = lif.euler(lif.V_reset)

    subplot(3,2,1), plot(time,potential), title('const h')
    subplot(3,2,2), plot(lif.h)

    lif.reset_spikes()
    lif.set_depolarizing_h();
    time, potential = lif.euler(lif.V_reset)

    subplot(3,2,3), plot(time,potential), title('depolarizing h')
    subplot(3,2,4), plot(lif.h)

    lif.reset_spikes()
    lif.set_hyperdepolarizing_h();
    time, potential = lif.euler(lif.V_reset)

    subplot(3,2,5), plot(time,potential), title('depolarizing h')
    subplot(3,2,6), plot(lif.h)


    #subplot(2,3,2), plot(lif.k), title('k')
    #subplot(2,3,3), plot(time,lif.i_stim), title('I_stim')
    #subplot(2,3,4), plot(lif.stim), title('stim')
    #subplot(2,3,5), plot(lif.h),
    show()


if __name__ == '__main__':
    #plot_fpt()
    #try_monte_carlo()
    #plot_monte_carlo_fpt()
    #try_pde()
    #plot_three_h()
    #compare_pde_mc_fpt()
    try_opt()
