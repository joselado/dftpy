from __future__ import print_function

import numpy as np
from . import calculations
from scipy.optimize import minimize

def volume(c0,n=0.3,relax=False,tol=0.005):
  """Optimize the volume of a certain calculation, returns 
  a new calculation with the optimized volume"""
  def f(x):
    """Function to minimize"""
    c = c0.copy() # copy calculation
    c.structure.expand(x[0]) # expand the structure
#    c.code = "QE" # use quantum Espresso
    e = c.gs_energy() # return the GS energy
    print("Computed energy for an expansion",x,e)
    return e # return the energy
  x = minimize(f,[1.0],bounds=[(1.0-n,1.0+n)],method='SLSQP',
                   options={"ftol":tol,"eps":tol/10.,}).x[0]
  print("Optimized",x)
  c = c0.copy() # create calculation
  c.structure.expand(x) # set this geometry
  return c
  



def lattice(c0,v=[2],n=0.3,relax=False,tol=0.001):
  """Optimize the lattice vectors of a certain calculation, returns 
  a new calculation with the optimized lattice"""
  nv = len(v) # number of variables
  def f(x):
    """Function to minimize"""
    c = c0.copy() # copy calculation
    dx = [1.,1.,1.] # initialize
    for i in range(nv): dx[v[i]] = x[i] # replace
    c.structure.expand(dx) # expand the structure
#    c.code = "QE" # use quantum Espresso
    e = c.gs_energy() # return the GS energy
    print("Computed energy for an expansion",dx,e)
    return e # return the energy
  bounds = [(1.0-n,1.0+n) for i in range(nv)] # initial bounds
  x0 = [1. for i in range(nv)] # initial guess
  x = minimize(f,x0,bounds=bounds,method='SLSQP',
                   options={"ftol":tol,"eps":tol/10.,}).x[0]
  print("Optimized",x)
  c = c0.copy() # create calculation
  c.structure.expand(x) # set this geometry
  return c



