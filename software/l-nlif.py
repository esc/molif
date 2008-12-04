from pylab import *
from scipy import *

# TODO:
# Add noise, at the moment it will have the same response
# the noise should be gaussian.
# what to do about the stimulus what should it be
# what would be a linear filter, a.k.a. kernel
# what about the I_hist, i don't even know what this should look like.


class lnlif:

    t_max = 10000

    # leak reversal potential : milli Volt
    V_leak = -80 
    # membrane leak conductance : siemens
    g =  1.e-4
    # reset voltage : milli Volt
    V_reset = -80
    # potential voltage : milli Volt
    V_threshold = -54

    # model of the stimulus as a vector : milli Ampere
    stim = zeros(t_max) 
    #noise = rand(t_max)

    # spationtemporal linear kernel
    # in this case modeled as a difference of gaussians
    n = stats.distributions.norm
    x = arange(-10,10,0.05)
    pos = stats.distributions.norm_gen.pdf(n,x,loc=0)
    neg = stats.distributions.norm_gen.pdf(n,x,loc=-4,scale=2)
    k = ( pos - neg)

    # now convolve the linear filter k with the stimulus
    # this really should be of length t_max
    # check that the convolution does indeed return this.
    # may need to pad out the linear kernel

    i_stim = convolve(stim,k, mode='same')
    spikes = []

    def set_rand_input(self):
        """ random input current """
        self.stim =rand(self.t_max)

    def set_const_input(self,current):
        """ constant input current """
        self.stim[:] = current

    def set_convolved_input(self):
        """ setup the input convolved with the linear filter"""
        self.i_stim = 0.2 + convolve(self.stim,self.k, mode='same')

    def integrate(self,x,t):
        return -self.g*(x - self.V_leak) + self.i_stim[t] #+ i_hist(t)
        
  
    def i_hist(self,t):
        # in the case where the temporal basis functions are not knowen
        # returning zero here brings us back to the more standard lif
        # model
        # todo this needs to be implemented
        return 0

    def euler(self, x_0, dt):
        """ euler method for for solving ODEs
        lif - a class that implements the integrate(x,t) function
        x_0 initial value for x
        t_max maximum time
        dt change in time """
        potential = zeros(self.t_max/dt)
        potential[0] = x_0 
        time = arange(0,self.t_max,dt)
        for i in xrange (1,int(self.t_max/dt)):
            if potential[i-1] >= self.V_threshold:
                self.spikes.append(i)
                potential[i] = self.V_reset;
            else:
                potential[i] = potential[i-1] + lif.integrate(potential[i-1],time[i-1]) * dt #+ noise[i]
        return (time, potential)

lif = lnlif()
lif.set_const_input(0.01);
lif.i_stim = lif.stim
# lif.set_convolved_input();


time, potential = lif.euler(-80,1)

print lif.spikes

subplot(2,3,1), plot(time,potential), title('output')
subplot(2,3,2), plot(lif.k), title('k')
subplot(2,3,3), plot(time,lif.i_stim), title('I_stim')
subplot(2,3,4), plot(lif.stim), title('stim')
#subplot(2,3,5), plot(lif.h),
show()


# notes:
# integrate for length of stimulus
