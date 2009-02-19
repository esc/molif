#!/usr/bin/env python
#coding=utf-8

# Copyright 2008 Valentin 'esc' Haenel
#
# This program is free software. It comes without any warranty, to
# the extent permitted by applicable law. You can redistribute it
# and/or modify it under the terms of the Do What The Fuck You Want
# To Public License, Version 2, as published by Sam Hocevar. See
# http://sam.zoy.org/wtfpl/COPYING for more details. 

from util import * 
from model import *
from numpy import nan, histogram, array , rot90, arange, isnan


def mc_P_vt_fpt():
    """ try to calculate both density and first passage time using
    monte carlo method """
    lif = lif_setup()
    lif.noise = True

    num_replications = 5000
    final_t_max = 400

    pots = zeros((num_replications,final_t_max))
    pots[:,:] = nan
    fpt = zeros(final_t_max)

    traces = []

    for i in xrange(num_replications):
        if i%100 == 0 : print i
        time, potential = \
        lif.euler(lif.V_reset,quit_after_first_spike=True)
        # now bin the potential 
        #for j in xrange(len(potential)):
        #    pots[i,j] = potential[j]

        pots[i,:len(potential)] = potential
        #and the fpt
        traces.append(potential)
        fpt[len(time)] += 1

    # make one histogram for each time step
    # WARNING THIS FITS EXACTLY THE PARAMETERS OF PDE SOLVER
    # DO NOT TOUCH THIS !!! EVER !!!
    V_range=arange(-3,1,4./500)
    P_vt = [ histogram(pots[~isnan(pots[:,i]),i],bins=V_range)[0] for i in \
            xrange(final_t_max)]
    P_vt = array(P_vt)
    # chop off top and bottom row, cause, they skew the colors
    # rot 90, to prepare for imshow()
    P_vt = rot90(P_vt)[1:-2,:]
    P_vt = (1.0* P_vt)/num_replications

    return P_vt, traces, fpt



@print_timing
def compute_mc_fpt(reps=500,t_max=400):
    """ compute only first passage time using monte carlo """
    print "computing monte carlo based first passage time now"
    print "using" , reps,"replications "
    lif = lif_setup()
    lif.noise = True

    n_steps = len(lif.stim)
    fpt = zeros(n_steps)
    spike_times = []

    for i in xrange(reps):
        if i%100 == 0 : print i
        time, potential = \
        lif.euler(lif.V_reset,quit_after_first_spike=True)
        if len(time) < n_steps:
            fpt[len(time)] += 1
            spike_times.append(time[-1])

    time = arange(n_steps)*lif.dt

    return time, fpt/reps/lif.dt, array(spike_times)


if __name__ == '__main__':
    pass
