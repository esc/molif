from pylab import *
from scipy import *

t_max = 10000

# leak reversal potential : volt
V_leak = -10
# membrane leak conductance
g = 0.01
# reset voltage, should be < 1
V_reset = -40
# threshold voltage
V_threshold = 15

# model of the stimulus as a vector
# TODO this should be fixed
stim = rand(t_max)
#noise = rand(t_max)

# spationtemporal linear kernel
# in this case modeled as a difference of gaussians
n = stats.distributions.norm
x = arange(-10,10,0.2)
pos = stats.distributions.norm_gen.pdf(n,x,loc=0)
neg = stats.distributions.norm_gen.pdf(n,x,loc=-5,scale=1.5)
k =  0.2 *( pos - neg)

# now convolve the linear filter k with the stimulus
# this really should be of length t_max
# check that the convolution does indeed return this.
# may need to pad out the linear kernel

I_stim = convolve(stim,k, mode='same')


def lnlif(x,t):
    ''' linear - noisy leaky integrate and fire model '''
    if x >= V_threshold:
        return V_reset
    else:
        return -g*(x - V_leak) + I_stim[t] #+ I_hist(t)
    
def I_hist(t):
    # in the case where the temporal basis functions are not knowen
    # returning zero here brings us back to the more standard LIF
    # model
    # TODO this needs to be implemented
    return 0


def euler(f_func, x_0, t_max, dt):
    """ euler method for for solving ODEs
    f_func - a method pointer
    x_0 initial value for x
    t_max maximum time
    dt change in time """
    result = zeros(t_max/dt)
    result[0] = x_0 
    time = arange(0,t_max,dt)
    for i in xrange (1,int(t_max/dt)):
        result[i] = result[i-1] + f_func(result[i-1],time[i-1]) * dt #+ noise[i]
    return (time, result)


time, result = euler(lnlif,-40,t_max,1)

subplot(2,2,1), plot(time,result), title('output')
subplot(2,2,2), plot(k), title('k')
subplot(2,2,3), plot(time,I_stim), title('I_stim')
subplot(2,2,4), plot(stim), title('stim')
show()

# notes:
# integrate for length of stimulus
