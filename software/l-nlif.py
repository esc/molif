#!/usr/bin/env python
#coding=utf-8

from pylab import *
from scipy import *
from scipy import stats, sparse, linsolve


class lnlif:

    def __init__(self, t_max=5000, dt=0.1):
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

    def pde_spike_train(self):
        """ compute the density evolution for all isis """

        # collection list for all densities
        densities = []

        for i in xrange(1,len(self.lif.spikes)):
            # for each ISI
            start = self.lif.spikes[i-1]
            end   = self.lif.spikes[i]
            if self.debug : print start
            if self.debug : print end

            density =  self.pde_interval(start,end)
            densities.append(density)
        
        return densities

    def pde_interval(self,start,end):
        """ compute the density evolution for the given interval """
        # length of the time interval
        U = end-start
        # final density matrix for this interval
        density = zeros((self.W,U))
        # set initial value
        density[self.V_reset_index,0] = 1
        for t in xrange(U-1):
            # V_rest as defined by lif
            V_rest = self.lif.V_rest(start+t);
            # A B and C diagonals of Lambda
            A = zeros(self.W-1)
            B = zeros(self.W)
            C = zeros(self.W-1)
            # zeroout the target
            #FIXME technically this could probably be skipped
            #since the values will be overwritten anyway
            self.beta[:] = 0
            # first we set the top boundary conditions
            # P(V_th,t) = 0 
            # i.e. the top row of the matrix
            B[0] = 1
            A[0] = 0
            self.beta[0] = 0
            # now we go and fill the matrix in
            a = self.a
            b = self.b
            c = self.c
            u = self.u
            w = self.w
            for j in xrange(self.W-2):
                b = self.lif.g * ((self.V_values[j]+self.V_values[j+1])/2. - V_rest)
                A[j+1] = -(2*a*u + b*w*u)
                B[j+1] = ((4*a*u) - (2*c*w**2*u) + (4*w**2))
                C[j] = -(2*a*u - b*w*u)
                self.beta[j+1] = (2*a*u + b*w*u) * density[j+2,t] + \
                        (-4*a*u + 2*c*w**2*u + 4*w**2) * density[j+1,t] + \
                        (2*a*u - b*w*u) * density[j,t]

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

            # FIXME this should be doable w/o for loop dumbass
            #for k in range(len(chi)):

            density[:,t+1] = chi[:]

            if self.debug: print "density: ", density

        return density

def try_pde():
    lif = lif_setup()

    pde = pde_solver(lif,500,-3.0,debug=False)

    #d = pde.pde_spike_train()
    #imshow(flipud(d[0]))
    d = pde.pde_interval(0,400)
    imshow(flipud(d))
    

    #for i in range(len(d)):
    #     subplot(4,len(d)/4,i), imshow(d[i])
    colorbar()

    P_vt = d

    figure()
    fpt = diff(P_vt.sum(axis=0))
    plot(fpt)
    show()


    return d

def compute_fpt():
    """ compute and the first passage time """
    print "computing partial differentail equation based first \
    passage time now"
    lif = lif_setup()
    pde = pde_solver(lif,500,-3.0,debug=False)
    P_vt = pde.pde_interval(0,400)
    fpt = abs(diff(P_vt.sum(axis=0)))
    
    return P_vt, fpt

def plot_fpt():
    P_vt, fpt = compute_fpt()
    imshow(flipud(P_vt))
    colorbar()
    figure()
    plot(fpt)
    show()

def lif_setup():

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
    lif.euler(lif.V_reset,quit_after_first_spike=True)

    return lif

def try_monte_carlo():

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

def monte_carlo_fpt(reps=500,t_max=400):
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

    mc_fpt = monte_carlo_fpt(reps=10000)
    P_vt, pde_fpt = compute_fpt()
    plot(mc_fpt,'r')
    plot(pde_fpt,'g')
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
    #plot_fpt()
    #try_monte_carlo()
    #plot_monte_carlo_fpt()
    #try_pde()
    #plot_three_h()


    compare_pde_mc_fpt()
