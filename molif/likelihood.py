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

def mle(variables,lif,V_lb,W,spikes,ids):
    """ maximum lkelihood function to give to optimizer """
    print "maximum likelihood function called"

    lif.g = variables[0]
    lif.V_leak = variables[1]
    lif.V_reset = variables[2]
    # ignoring k and h for now

    lif.spikes = spikes

    # what about dt?

    p = pde_solver(lif,W,V_lb)
    likelihood = p.compute_product_fpt()
    del p
    return likelihood

@print_timing
def try_opt():
    print "trying optimizer"
    # create the neuron
    lif = lif_setup()
    # generate some spikes
    time, potential = lif.euler(lif.V_reset)

    # make the variables vector
    # [g,V_leak,V_reset,k,h]
    variables = zeros(3+len(lif.k)+len(lif.h))
    id_g = 0
    id_V_leak = 1
    id_V_reset = 2
    #id_k = 3
    #id_h = 3 + len(k)
    ids = (id_g,id_V_leak,id_V_reset)
    variables[id_g] = lif.g
    variables[id_V_leak] = lif.V_leak
    variables[id_V_reset] = lif.V_reset
    #variables[id_k:id_h] = lif.k
    #variables[id_h:-1] = lif.h

    # make the fixed tuple 
    fixed = (lif,-3.0,500,lif.spikes,ids)

    print mle(variables,lif,-2.0,100,lif.spikes,ids)

    #xopt = optimize.fmin(mle,variables,fixed)

    #print xopt


