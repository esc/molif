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

    # leak reversal potential : volt
    V_leak = -80
    # membrane leak conductance
    g = 0.00001
    # reset voltage, should be < 1
    V_reset = -80
    # threshold voltage
    V_threshold = -54

    # model of the stimulus as a vector
    # TODO this should be fixed
    stim = rand(t_max)
    #noise = rand(t_max)

    # spationtemporal linear kernel
    # in this case modeled as a difference of gaussians
    n = stats.distributions.norm
    x = arange(-10,10,0.05)
    pos = stats.distributions.norm_gen.pdf(n,x,loc=0)
    neg = stats.distributions.norm_gen.pdf(n,x,loc=-4,scale=2)
    k = 0.7 * ( pos - neg)

    # now convolve the linear filter k with the stimulus
    # this really should be of length t_max
    # check that the convolution does indeed return this.
    # may need to pad out the linear kernel

    I_stim = 0.2 + convolve(stim,k, mode='same')
    spikes = []


    def integrate(self,x,t):
        if x >= self.V_threshold:
            self.spikes.append(t)
            return self.V_reset
        else:
            return -self.g*(x - self.V_leak) + self.I_stim[t] #+ I_hist(t)
        
  
    def I_hist(self,t):
        # in the case where the temporal basis functions are not knowen
        # returning zero here brings us back to the more standard LIF
        # model
        # TODO this needs to be implemented
        return 0


def euler(lif, x_0, t_max, dt):
    """ euler method for for solving ODEs
    lif - a class that implements the integrate(x,t) function
    x_0 initial value for x
    t_max maximum time
    dt change in time """
    result = zeros(t_max/dt)
    result[0] = x_0 
    time = arange(0,t_max,dt)
    for i in xrange (1,int(t_max/dt)):
        result[i] = result[i-1] + lif.integrate(result[i-1],time[i-1]) * dt #+ noise[i]
    return (time, result)

lif = lnlif()

time, result = euler(lif,-80,lif.t_max,1)

print lif.spikes

subplot(2,2,1), plot(time,result), title('output')
subplot(2,2,2), plot(lif.k), title('k')
subplot(2,2,3), plot(time,lif.I_stim), title('I_stim')
subplot(2,2,4), plot(lif.stim), title('stim')
show()


# notes:
# integrate for length of stimulus
