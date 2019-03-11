from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

import structure
import matplotlib.pyplot as plt

def get(q):
  def get_struct(s=2.7):
    """Get a structure for Fe"""
    a123 = np.array([[1,1,-1],[1,-1,1],[-1,1,1]])*s # lattice vectors
    atoms = [["Fe",[0.,0.,0.]]] # first atom
    st = structure.create_structure(a123,atoms) # create structure
    return st
  
  st = get_struct() # Fe structure
  
  import calculations # calculations module
  c = calculations.Calculation(st,spin=True,m=[0.,0.,-1.]) 
  c.kmesh = [1,1,1]
  c.set_spiral(q)
  e0 = c.gs_energy() # energy
  return e0

f = open("SWEEP.OUT","w")

for x in np.linspace(0.,2.,10):
  for y in np.linspace(0.,2.,10):
    print("Doing",x,y)
    e = get([x,y,0])
    f.write(str(x)+"  "+str(y)+"  "+str(e)+"\n")
    f.flush()


f.close()
