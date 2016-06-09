# -*- coding: utf-8 -*-
"""
Created on Sun May 22 12:15:34 2016

@author: S.Wystub
"""

import matplotlib.pyplot as plt
import numpy as np

#This code is used to create the compactness over thickness plot for a anisotropic gravastars, see this paper for more: https://arxiv.org/abs/0706.1513
#By commenting out lines 84-86 and un-commenting line 92, the pogram can be used to derive the other plots for anisotropic gravastars in the paper
#Oh god why


def compute_density(r,r1,r2,density0):
  ''' compute the equation of state, i.e. eq 23 '''
  if r <= r1 and r >= 0:
    return density0
  elif r1 >= r or r2 <= r:
    return 0

  #defining density function
  a = (2*density0)
  b = (-3*density0)*(r2 + r1)
  c = 6*density0*r1*r2
  d = density0*(r2**3 - 3*r1*r2**2)

  return (a*r**3 + b*r**2 + c*r + d)/(r2 - r1)**3

#defining core program: parse r1, r2, m(r2) and a stepsize for evaluation, return mass, pressures, metric functions
def get_grr(r1, r2, M, stepsize): 
  ''' compute 1-2m/r where m is gravitational mass'''
  #creating an array to evaluate the functions at every point
  r = np.arange(0, 2*r2, stepsize)
  #avoid origin
  r[0]=1.e-16
  
  #fixing rho_0 in terms of M
  density0 = M*((r2 - r1)**3/(4*np.pi))*np.power(((r2**6 - r1**6)/3) - (3*(r2 + r1)*(r2**5 - r1**5)/5) + (3*r1*r2*(r2**4 - r1**4)/2) + ((r2**3 - 3*r2**2*r1 )*(r2**3 - r1**3)/3) + ((r1**3*(r2 - r1)**3)/3), -1) 

  vfunc = np.vectorize(compute_density)
  dens=vfunc(r,r1,r2,density0) 
  mass=np.cumsum(dens*r**2)*np.pi*4*stepsize 
    
  return (1-2*mass/r)

def is_neg(grr):
  ''' return if grr<0 '''
  is_zero=len(np.where(grr<=0)[0])
  if is_zero>0:
    return 1
  else:
    return 0
  
if __name__ == '__main__':

  #parsing arguments (rescaled since the bash script sequence only allows integers, e.g. parse 100 if you want a mass of 1)
#  r1 = float(sys.argv[1])
#  r2 = float(sys.argv[2])
#  M = float(sys.argv[3])
#  stepsize = float(sys.argv[4])
  
  stepsize=0.1
  r2=2.0 
  r1_max=r2-0.01
  m_max=3.0
  mu=[]
  delta=[]
  for r1 in np.linspace(0.0,r1_max,21):
    #print r1
    for M in np.linspace(0.01,m_max,301): 
      grr=get_grr(r1, r2, M, stepsize)
      if not is_neg(grr):
        mu.append(M/r2)
        delta.append((r2-r1)/M)

plt.scatter(delta,mu)
plt.xlim([0.0,20])
plt.ylim([0,0.5])
plt.show()
