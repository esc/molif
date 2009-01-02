#!/usr/bin/env python
#coding=utf-8

from pylab import *
from scipy import *

class lnlif:

    def __init__(self):
        # a place to store the spikes
        self.spikes = [0]

        # number of time slots to integrate until
        self.t_max = 1000
        self.dt = 1

        # leak reversal potential : milli Volt
        self.V_leak = 0.9
        # membrane leak conductance : siemens
        self.g = 1.e-2
        # reset voltage : milli Volt
        self.V_reset = 0
        # potential voltage : milli Volt
        self.V_threshold = 1

        # model of the stimulus as a vector : milli Ampere
        self.stim = zeros(self.t_max) 

        # st.d. of the gaussian white noise process.
        self.sigma = 0.02
        self.noise = False
        self.h_scale = 0.5

        # spationtemporal linear kernel
        # in this case modeled as a difference of gaussians
        n = stats.distributions.norm
        x = arange(-10,10,0.05)
        pos = stats.distributions.norm_gen.pdf(n,x,loc=0)
        neg = stats.distributions.norm_gen.pdf(n,x,loc=-4,scale=2)
        self.k = ( pos - neg)
        

    def set_rand_input(self):
        """ random input current """
        self.stim =rand(self.t_max)

    def set_const_input(self,current):
        """ constant input current """
        self.stim[:] = current

    def set_depolarizing_h(self):
        self.h = self.h_scale * 1/exp(self.get_time())

    def set_const_h(self):
        self.h = zeros(len(self.get_time()))

    def set_hyperdepolarizing_h(self):
        self.h =  self.h_scale * -1/exp(self.get_time())

    def get_time(self):
        return arange(0,self.t_max/self.dt,self.dt)

    def reset_spikes(self):
        self.spikes = [0]


    def set_convolved_input(self):
        """ setup the input convolved with the linear filter"""
        # now convolve the linear filter k with the stimulus
        # this really should be of length t_max
        # check that the convolution does indeed return this.
        # may need to pad out the linear kernel
        self.i_stim = 0.2 + convolve(self.stim,self.k, mode='same')

    def integrate(self,x,t):
        return -self.g*(x - self.V_leak) + self.i_stim[t] + self.i_hist(t)
        
  
    def i_hist(self,t):
        # in the case where the temporal basis functions are not knowen
        # returning zero here brings us back to the more standard lif
        # model
        return self.h[t-self.spikes[-1]]

    def add_noise(self, dt):
        """ gaussian additive white noise"""
        if (self.noise):
            return self.sigma * sqrt(dt) * random.randn();
        else:
            return 0


    def euler(self, x_0):
        """ euler method for for solving ODEs
        lif - a class that implements the integrate(x,t) function
        x_0 initial value for x
        t_max maximum time
        dt change in time """
        self.potential = zeros(self.t_max/self.dt)
        self.potential[0] = x_0 
        self.time = arange(0,self.t_max,self.dt)
        for i in xrange (1,int(self.t_max/self.dt)):
            if self.potential[i-1] >= self.V_threshold:
                self.spikes.append(i)
                self.potential[i] = self.V_reset;
            else:
                self.potential[i] = self.potential[i-1] + \
                self.integrate(self.potential[i-1],self.time[i-1]) * self.dt + \
                self.add_noise(self.dt)
        return (self.time, self.potential)

    def V_rest(self,t):
        return self.V_leak + 

def pde_solver(lif,W,U,V_lb):
    # lif provides voltage trace, time, stimulus, threshold, an
    # estimate of the spatio temporal linear kernel and an estimation
    # of the h function. 
    # W is the number of points at which V will be calculated
    # U is the number of 
    # V_lb is the lower bound on the voltage discretization, the V_th
    # from lif is the uppper bound.

    # current p is a vector of length W
    # index 0 is V_lb and index end is V_th
    # the initial value is zero everywhere except V_reset

    # here we just preallocate
    # this should be a matrix
    current_p = zeros(W) 
    # these are the values 
    current_p_values = zeros(W)

    # Allocate a sparse matrix for later use
    A = scipy.sparse.lil_matrix()
    # set the variables needed for filling in the matrix
    # these will not change 
    a = 1/2
    c = lif.g
    # this one will change for each reinitialization of the 
    # matrix A, even for each entry.
    b = 0

    
    # for each isi
        # for each time step
            # calculate V_rest
            # initialize A
            # initialize b
            # solve Ax = b
            # store x in current_p matrix
   
    

    # extract next interspike interval
    
    pass



def plot_three_h():

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
    plot_three_h()
