from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

import structure


def get_struct(s=9.3/2.0):
  """Get a structure for silicon"""
  a123 = np.array([[1,1,0],[1,0,1],[0,1,1]])*s # silicon lattice vectors
  atoms = [["Pb",[0.,0.,0.]]] # first atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # silicon structure

import calculations # calculations module

# quantum espresso
c1 = calculations.Calculation(st,code="Elk",soc=False) # create a DFT calculation
ks1,es1 = c1.bandstructure() # get bandstructure
ks1 = ks1/max(ks1)


# Elk
c2 = calculations.Calculation(st,code="Elk",soc=True) # create a DFT calculation
ks2,es2 = c2.bandstructure() # get bandstructure
ks2 = ks2/max(ks2)


import matplotlib.pyplot as plt

plt.scatter(ks1,es1,c="red",label="No SOC")
plt.scatter(ks2,es2,c="blue",label="SOC")

plt.legend()
plt.show()

