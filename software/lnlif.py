#!/usr/bin/env python
#coding=utf-8

from pylab import imshow, figure, plot, show
from numpy import arange, zeros, sqrt, random, diff, cumsum
from scipy import stats, sparse, linsolve, optimize
import time
import cProfile


def print_timing(func):
    """ decorator to estimate the timing of methods """
    def wrapper(*arg, **kwargs):
        t1 = time.time()
        res = func(*arg , **kwargs)
        t2 = time.time()
        print '%s took %0.3f ms' % (func.func_name, (t2-t1)*1000.0)
        return res
    return wrapper


#TODO pylint it
#     the definition of I_hist isn't fully implemented
#     fix the implementation of the kernel


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

class pde_solver():
    def __init__(self,lif,W,V_lb,debug=False):
        """ compute the density evolution for a lif model

            Arguments:
            lif  - and instance of lnlif
            W    - the number of intervals to split lif.t_max into
            V_lb - the lower bound on the voltage discretization
            debug- when True will output lots of debugging

            lif provides most of the variables for this function.

        """

        self.lif = lif
        self.W = W
        self.V_lb = V_lb
        self.debug = debug
        # Lambda is our tridiagonal matrix in each iteration step
        self.Lambda = sparse.lil_matrix((self.W,self.W))
        # coefficients of the terms in the Fokker-Planck equation
        self.a = self.lif.sigma**2/2
        self.c = self.lif.g
        self.b = 0
        # length of time interval in discretization
        self.u = self.lif.dt
        # voltage upper bound
        self.V_max = self.lif.V_threshold;
        # voltage lower bound
        self.V_min = self.V_lb
        # voltage range 
        self.V_range = self.V_max - self.V_min
        # lenth of a voltage interval
        self.w = self.V_range/self.W
        # mapping of possible values for voltage note: len(V_values) = W
        self.V_values = arange(self.V_min,self.V_max,self.w)
        if self.debug: print "V_values: " , self.V_values
        # this is the index of the value thats closest to V_reset
        self.V_reset_index = abs(self.V_values - lif.V_reset).argmin()  
        if self.debug: print "V_reset_index" , self.V_reset_index
        # Lambda * chi = beta
        self.beta = zeros(W)

    def compute_product_fpt(self):
        """ compute the product of all fpts for all spike intervals
        this is the maximum likelihood
        """
        print "computing fpt for all intervals, total: ", \
        len(self.lif.spikes)
        likelihood = 0
        for i in xrange(1,len(self.lif.spikes)):
            if i%10==0 : print "interval: " , i
            # for each ISI
            start = self.lif.spikes[i-1]
            end   = self.lif.spikes[i]
            if self.debug : 
                print "start: ", start
                print "end:" , end
                print "spike_times", self.lif.spikes

            density =  self.pde_interval(start,end)
            likelihood += self.density_to_fpt(density)[-1]
            del density
        
        return likelihood

    def density_to_fpt(self,density):
        """ turn the density into an fpt
        Differentiat w.r.t. time the integral w.r.t. Potential

        """
        return diff(density.sum(axis=0)) * -1.0

    #@print_timing
    def pde_interval(self,start,end):
        """ compute the density evolution for the given interval """
        # length of the time interval
        U = end-start
        # final density matrix for this interval
        density = zeros((self.W,U))
        # set initial value
        if self.debug : 
            print "W" , self.W
            print "U" , U
            print density 
        density[self.V_reset_index,0] = 1

        for t in xrange(U-1):
            # V_rest as defined by lif
            V_rest = self.lif.V_rest(start+t);
            # A B and C diagonals of Lambda
            A = zeros(self.W-1)
            B = zeros(self.W)
            C = zeros(self.W-1)
            # first we set the top boundary conditions
            # P(V_th,t) = 0 
            # i.e. the top row of the matrix
            B[0] = 1
            A[0] = 0
            self.beta[0] = 0
            # now we go and fill the matrix in
            a,b,c,u,w = self.a, self.b, self.c, self.u, self.w

            b_array = self.lif.g * ((self.V_values[:-2] + \
                self.V_values[1:-1])/2 -V_rest)
            
            #here we scrapped the for loop, yes!
            #but scrapped readability, no!
            A[1:]= -(2*a*u + b_array*w*u)

            B[1:-1] = ((4*a*u) - (2*c*w**2*u) + (4*w**2))

            C[:-1]= -(2*a*u - b_array*w*u)
            
            self.beta[1:-1] = (2*a*u + b_array[:]*w*u) * density[2:,t] + \
                (-4*a*u + 2*c*w**2*u + 4*w**2) * density[1:-1,t] +   \
                (2*a*u - b_array[:]*w*u) * density[:-2,t]


            # now we need to fill in the last row, which are the lower
            # boundary conditions, i.e. P(V_lb,t) = 0
            C[self.W-2] = 0
            B[self.W-1] = 1
            self.beta[self.W-1] = 0
            # now we set the diagonals of the tridiagonal matrix
            self.Lambda.setdiag(A,1)
            self.Lambda.setdiag(B)
            self.Lambda.setdiag(C,-1)

            
            if self.debug: 
                print "A :" , A
                print "B :" , B
                print "C :" , C
                print "Lambda: " , self.Lambda.todense()
                print "beta: " , self.beta
            
            chi = linsolve.spsolve(self.Lambda,self.beta)

            if self.debug: 
                print "chi:" , chi
                print "sumchi", chi.sum()

            density[:,t+1] = chi[:]

            if self.debug: print "density: ", density

        return density

def try_pde():
    """ try pde solver without simulating a single spike of the neuron
    """
    lif = lif_setup()
    pde = pde_solver(lif,500,-3.0,debug=False)
    d = pde.pde_interval(0,400)

    imshow(flipud(d))
    colorbar()
    P_vt = d
    figure()
    fpt = diff(P_vt.sum(axis=0))
    plot(fpt)
    show()

    return d

@print_timing
def compute_pde_fpt():
    """ compute the first passage time using pde"""
    print "computing partial differentail equation based first \
    passage time now"
    lif = lif_setup()
    p = pde_solver(lif,500,-3.0,debug=False)
    density = p.pde_interval(0,400)
    fpt = p.density_to_fpt(density)
    return density, fpt

def plot_fpt():
    """ plot density evolution and and first passage time. """
    P_vt, fpt = compute_pde_fpt()
    imshow(flipud(P_vt))
    colorbar()
    figure()
    plot(fpt)
    show()

def lif_setup():
    """ initialize nueron model with default parameters """
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

def plot_monte_carlo_fpt():
    mc_fpt = monte_carlo_fpt(reps=10000)
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