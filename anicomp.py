# -*- coding: utf-8 -*-
"""
Created on Sun May 22 12:15:34 2016

@author: S.Wystub
"""

import scipy.integrate as integrate
import numpy as np
import copy
from scipy.misc import derivative
import sys

#This code is used to create the compactness over thickness plot for a anisotropic gravastars, see this paper for more: https://arxiv.org/abs/0706.1513
#By commenting out lines 84-86 and un-commenting line 92, the pogram can be used to derive the other plots for anisotropic gravastars in the paper
#Oh god why


#ignoring div0 errors
np.seterr(divide='ignore', invalid='ignore')

#parsing arguments (rescaled since the bash script sequence only allows integers, e.g. parse 100 if you want a mass of 1)
r1 = 0.1*float(sys.argv[1])
r2 = 0.1*float(sys.argv[2])
M = 0.01*float(sys.argv[3])
stepsize = float(sys.argv[4])

#defining core program: parse r1, r2, m(r2) and a stepsize for evaluation, return mass, pressures, metric functions
def core(r1, r2, M, stepsize):

	#creating an array to evaluate the functions at every point
	inputarray = np.arange(0, 2*r2, stepsize)
	
	#fixing rho_0 in terms of M
	density0 = M*((r2 - r1)**3/(4*np.pi))*np.power(((r2**6 - r1**6)/3) - (3*(r2 + r1)*(r2**5 - r1**5)/5) + (3*r1*r2*(r2**4 - r1**4)/2) + ((r2**3 - 3*r2**2*r1 )*(r2**3 - r1**3)/3) + ((r1**3*(r2 - r1)**3)/3), -1)
	
	#defining density function
	a = (2*density0)/(r2 - r1)**3
	b = (-3*density0)*(r2 + r1)/(r2 - r1)**3
	c = 6*density0*r1*r2/(r2 - r1)**3
	d = density0*(r2**3 - 3*r1*r2**2)/(r2 - r1)**3
	density = lambda r: np.piecewise(r, [(r<r1), (r>= r1) & (r<=r2), (r>r2)], [lambda r: density0, lambda r: a*r**3 + b*r**2 + c*r + d, lambda r: 0])
	
	#defining mass function
	def mass(r):
		dummyarray = copy.copy(r)
		for x in np.nditer(dummyarray, op_flags=['readwrite']):
			x[...] = integrate.quad(lambda t: 4*np.pi*density(np.array([t]))*t**2, 0, x)[0]
		return dummyarray
		
	#fixing alpha (derived by dp_r/dr=1 at maximum)
	alpha = 2.2135
	
	#defining radial pressure
	def pr(r):
		return (density(r)**2/density0)*(alpha - (1 + alpha)*(density(r)/density0)**2)
	
	#defining pr'(r)
	def dpr(r):
		return derivative(pr, r, dx=1e-6)
	
	#defining tangential pressure
	def pt(r):
		return ((r*dpr(r)/2) + pr(r) + ((density(r) + pr(r))/2)*((4*np.pi*pr(r)*r**3 + mass(r))/(r-2*mass(r))))
	
	#defining metric function g_rr
	def grr(r):
		return (1-(2*mass(r))/r)
		
	#defining gamma(r) and g_tt, commented out for performance. Double integrals are tricky.
    #def gamma(r):
    #    dummyarray = copy.copy(r)
    #    for x in np.nditer(dummyarray, op_flags=['readwrite']):
    #        x[...] = integrate.quad(lambda t: (2*mass(np.array([t]))+8*np.pi*pr(np.array([t]))*t**3)/(t**2 -2*t*mass(np.array([t]))), 0, x)[0]
    #    return dummyarray
    
    #def gtt(r):
    #    return -(1-2*M/r2)*np.exp(gamma(r)-gamma(np.array([r2])))
    
    #not returning g_tt for performance reasons
	
	#checking g_rr for singularities (negative values), returning thickness/M, compactness for the plot
	#singular solutions are sorted out by returning compactness -1
	if np.nanmin(grr(inputarray)) > 0:
		return ((r2-r1)/M, M/r2)
	else:
		return ((r2-r1)/M, -1)
		

    
    #returning the functions characterizing the gravastar
	#not used here
    #return np.array([inputarray, density(inputarray), mass(inputarray), pr(inputarray), pt(inputarray), grr(inputarray)])
    
 
print(core(r1, r2, M, stepsize))
