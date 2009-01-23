#!/usr/bin/env python
#coding=utf-8

# Copyright 2008 Valentin 'esc' Haenel
#
# This program is free software. It comes without any warranty, to
# the extent permitted by applicable law. You can redistribute it
# and/or modify it under the terms of the Do What The Fuck You Want
# To Public License, Version 2, as published by Sam Hocevar. See
# http://sam.zoy.org/wtfpl/COPYING for more details. 

from numpy import arange, zeros, sqrt, random
from scipy import stats


class lnlif:

    def __init__(self, t_max=100, dt=0.1):
        """ initialize the L-NLIF model using max_time and dt
        DO NOT CHANGE t_max or dt AFTER INITIALIZATION """
        # number of time slots to integrate until
        self.t_max = t_max
        # discretization in time
        self.dt = dt
        # time vector
        self.time = arange(0,self.t_max/self.dt,self.dt)

        # a place to store the spikes
        self.spikes = [0]

        # leak reversal potential : milli Volt
        self.V_leak = 0.9
        # membrane leak conductance : siemens
        self.g = 1.e-2
        # reset voltage : milli Volt
        self.V_reset = 0
        # potential voltage : milli Volt
        self.V_threshold = 1
        # allocate for stimulus : milli Ampere
        self.stim = zeros(self.t_max/self.dt) 
        # st.d. of the gaussian white noise process.
        self.sigma = 0.02
        self.noise = False
        # scaling of the afterspike current waveform
        self.h_scale = 0.5

        #FIXME need some way of setting this
        # spationtemporal linear kernel
        # in this case modeled as a difference of gaussians
        # for starters just use a delta pulse, as the most simple
        # kernel, or even a short filter(couple of values)
        n = stats.distributions.norm
        x = arange(-10,10,0.05)
        pos = stats.distributions.norm_gen.pdf(n,x,loc=0)
        neg = stats.distributions.norm_gen.pdf(n,x,loc=-4,scale=2)
        self.k = ( pos - neg)
        

    def set_rand_input(self):
        self.stim =rand(self.t_max) * 0.5

    def set_const_input(self,current):
        self.stim[:] = current

    def set_depolarizing_h(self):
        self.h = self.h_scale * 1/exp(self.time())

    def set_const_h(self):
        self.h = zeros(len(self.time))

    def set_hyperdepolarizing_h(self):
        self.h =  self.h_scale * -1/exp(self.time())

    def get_time_vector(self):
        """ returns the tim indices """

    def reset_spikes(self):
        self.spikes = [0]

    def set_convolved_input(self):
        """ setup the input convolved with the linear filter"""
        # now convolve the linear filter k with the stimulus
        # FIXME this really should be of length t_max
        # check that the convolution does indeed return this.
        # may need to pad out the linear kernel
        self.i_stim = 0.2 + convolve(self.stim,self.k, mode='same')

    def integrate(self,x,t):
        """ definition of dV/dt """
        return -self.g*(x - self.V_leak) + self.i_stim[t] + self.i_hist(t)
  
    def i_hist(self,t):
        # in the case where the temporal basis functions are not knowen
        # returning zero here brings us back to the more standard lif
        # model
        return self.h[t-self.spikes[-1]]

    def add_noise(self):
        """ gaussian additive white noise"""
        if (self.noise):
            return self.sigma * sqrt(self.dt) * random.randn();
        else:
            return 0


    def euler(self, x_0, quit_after_first_spike=False):
        """ euler method for for solving ODEs
        lif - a class that implements the integrate(x,t) function
        x_0 initial value for x
        t_max maximum time
        dt change in time """
        #TODO needs also a t_max
        # i.e. a time to integrate until
        self.potential = zeros(self.t_max/self.dt)
        self.potential[0] = x_0
        self.reset_spikes()
        for i in xrange (1,int(self.t_max/self.dt)):
            if self.potential[i-1] >= self.V_threshold:
                self.spikes.append(i)
                self.potential[i] = self.V_reset;
                if quit_after_first_spike:
                    return (self.time[:i] , self.potential[:i])
            else:
                self.potential[i] = self.potential[i-1] + \
                self.integrate(self.potential[i-1],self.time[i-1]) * self.dt + \
                self.add_noise()
        return (self.time, self.potential)

    def V_rest(self,t):
        """ stationary point of the noiseless subthreshold dynamics
        """
        #return 0
        return self.V_leak + 1/self.g * (self.i_stim[t] + self.i_hist(t));

def lif_setup():
    """ initialize nueron model with dsome efault parameters """
    lif = lnlif(dt=0.01) # init model
    lif.set_const_input(0.5); # set constant input
    #lif.set_rand_input()
    lif.i_stim = lif.stim # setup stimulus
    # lif.set_convolved_input();
    lif.noise = False
    lif.set_const_h()
    lif.V_leak = 0.0
    #lif.g = 0
    lif.sigma = 0.1

    time, potential = \
    lif.euler(lif.V_reset,quit_after_first_spike=False)

    return lif


if __name__ == '__main__':
    print  "loaded molif.model, will simulate a LNLIF now "
    l = lif_setup()
    print l.__dict__

