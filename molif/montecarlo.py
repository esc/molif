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


def try_monte_carlo():
    """ try to calculate both density and first passage time using
    monte carlo method """
    lif = lif_setup()
    lif.noise = True

    num_replications = 500
    final_t_max = 400

    potentials = zeros((num_replications,lif.t_max))
    potentials[:,:] = NaN
    first_passage_time = zeros(lif.t_max)
    ioff()
    figure()
    hold(True)

    V_max = 1
    V_min = 0
    V_step = 0.005
    V_range=arange(V_min,V_max,V_step)

    for i in xrange(num_replications):
        print i
        time, potential = \
        lif.euler(lif.V_reset,quit_after_first_spike=True)
        plot(time,potential)
        # now bin the potential 
        for j in xrange(len(potential)):
            potentials[i,j] = potential[j]
        #and the fpt
        first_passage_time[len(time)] += 1
        if len(time) > final_t_max:
            final_t_max = len(time)

    print "final_t_max: " , final_t_max

    figure()

    # make one histogram for each time step
    y = [ histogram(potentials[:,i],bins=V_range)[0] for i in \
            xrange(final_t_max)]
    y = array(y)
   
    # chop off top and bottom row, cause, they skew the colors
    print shape(y)
    imshow(rot90(y)[1:-2,:])

    colorbar()
    figure()
    plot(first_passage_time[:final_t_max])
    show()

@print_timing
def monte_carlo_fpt(reps=500,t_max=400):
    """ compute only first passage time using monte carlo """
    print "computing monte carlo based first passage time now"
    print "using" , reps,"replications "
    lif = lif_setup()
    lif.noise = True

    fpt = zeros(t_max)

    for i in xrange(reps):
        if i%100 == 0 : print i
        time, potential = \
        lif.euler(lif.V_reset,quit_after_first_spike=True)
        fpt[len(time)] += 1

    return fpt/reps


