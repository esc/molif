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

from numpy import diff
from scipy import sparse, linsolve

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

            P_vt =  self.pde_interval(start,end)
            likelihood += self.P_vt_to_fpt(P_vt)[-1]
            del P_vt
        
        return likelihood

    def P_vt_to_fpt(self,P_vt):
        """ turn the density evolution into an first passage time
        Differentiat w.r.t. time the integral w.r.t. Potential

        """
        return diff(P_vt.sum(axis=0)) * -1.0

    #@print_timing
    def pde_interval(self,start,end):
        """ compute the density evolution for the given interval """
        # length of the time interval
        U = end-start
        # final density matrix for this interval
        P_vt = zeros((self.W,U))
        # set initial value
        if self.debug : 
            print "W" , self.W
            print "U" , U
            print P_vt
        P_vt[self.V_reset_index,0] = 1

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
            
            self.beta[1:-1] = (2*a*u + b_array[:]*w*u) * P_vt[2:,t] + \
                (-4*a*u + 2*c*w**2*u + 4*w**2) * P_vt[1:-1,t] +   \
                (2*a*u - b_array[:]*w*u) * P_vt[:-2,t]


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

            P_vt[:,t+1] = chi[:]

            if self.debug: print "P_vt: ", P_vt

        return P_vt

@print_timing
def compute_pde_fpt():
    """ compute the first passage time using pde"""
    print "computing partial differentail equation based first \
    passage time now"
    lif = lif_setup()
    p = pde_solver(lif,500,-3.0,debug=False)
    P_vt = p.pde_interval(0,400)
    fpt = p.P_vt_to_fpt(P_vt)
    return P_vt, fpt


if __name__ == '__main__':
    pass

