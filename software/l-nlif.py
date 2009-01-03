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

    def add_noise(self):
        """ gaussian additive white noise"""
        if (self.noise):
            return self.sigma * sqrt(self.dt) * random.randn();
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
                self.add_noise()
        return (self.time, self.potential)

    def V_rest(self,t):
        return self.V_leak + 1/self.g * (self.i_stim[t] + self.i_hist(t));

def pde_solver(lif,W,V_lb):
    # lif provides voltage trace, time, stimulus, threshold, an
    # estimate of the spatio temporal linear kernel and an estimation
    # of the h function. 
    # W is the number of points at which P(V,t) will be calculated
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
    Lambda = sparse.lil_matrix((W,W))
    # set the variables needed for filling in the matrix
    # these will not change 
    a = 1/2
    c = lif.g
    # this one will change for each reinitialization of the 
    # matrix A, even for each entry.
    b = 0
    # take the length of a time interval to be the 
    # dt of the neuron we used as input
    u = lif.dt


    V_max = lif.V_threshold;
    V_min = V_lb
    # now divide the difference into W equally spaced intervals
    V_range = V_max - V_min
    w = V_range/W
    V_values = arange(V_min,V_max,w)
    print "V_values: " , V_values
    # note: len(V_values) = W
    # this is the index of the value V_reset
    V_reset_index = abs(V_values - lif.V_reset).argmin()  
    print "V_reset_index" , V_reset_index

    # collect all densities for all ISIs
    densities = []
    # this will be the target for Ax = b
    beta = zeros(W)

    for i in xrange(1,len(lif.spikes)):
        # for each ISI
        start = lif.spikes[i-1]
        end   = lif.spikes[i]
        # this is our time discretization
        # because we may not have the value of the stimulus at
        # alternative values this may be different for each ISI
        # if there is noise
        U = end-start
        # note: U is also the number of rows in the final
        # density matrix to store the result of the density evolution.
        density = zeros((W,U))
        #here we set the left hand side boundary condition
        density[V_reset_index,0] = 1
        for t in xrange(U-1):
            # for each time step
            rest = lif.V_rest(start+t);
            # now we really need to initialize the target matrix
            # we do this stupidly first and then optimize
            # now we need the A B and C diagonals
            # A is the top diagonal
            A = zeros(W-1)
            # B is THE diagonal
            B = zeros(W)
            # C is the lower diagonal
            C = zeros(W-1)
            # zeroout the target
            #FIXME technically this could probably be skipped
            #since the values will be overwritten anyway
            beta[:] = 0
            # first we set the top boundary conditions
            # P(V_th,t) = 0 
            # i.e. the top row of the matrix
            B[0] = 1
            A[0] = 0
            beta[0] = 0
            # now we go and fill the matrix in
            for j in xrange(W-2):
                b = lif.g * (V_values[-j-1] - rest)
                A[j+1] = -(2*a*u + b*w*u)
                B[j+1] = ((4*a*u) - (2*c*(w**2)*u) + (4*w**2))
                C[j] = -(2*a*u + b - w * u)
                beta[j+1] = (2*a*u + b*w*u) * density[j,t] + \
                        (-4*a*u + 2*c*w**2*u + 4*w**2) * density[j+1,t] + \
                        (2*a*u - b*w*u) * density[j+2,t]

            # now we need to fill in the last row
            C[W-2] = 0
            B[W-1] = 1
            beta[W-1] = 0
            # now we setup the tridiagonal matrix
            Lambda.setdiag(A,1)
            Lambda.setdiag(B)
            Lambda.setdiag(C,-1)

            print "A :" , A
            print "B :" , B
            print "C :" , C
         
            print "Lambda: " , Lambda
            print "beta: " , beta
        
            chi = linsolve.spsolve(Lambda,beta)
            print chi

            # FIXME this should be doable w/o for loop dumbass
            for k in range(len(chi)):
                density[k,t+1] = chi[k]
            
        densities.append(density)

        
    return densities

def try_pde():

    lif = lnlif() # init model
    lif.set_const_input(0.01); # set constant input
    lif.i_stim = lif.stim # setup stimulus
    # lif.set_convolved_input();
    #lif.noise = True
    lif.set_const_h();

    time, potential = lif.euler(lif.V_reset)

    d = pde_solver(lif,100,-0.2)

    imshow(d[0])
    show()


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
    try_pde()
    #plot_three_h()
