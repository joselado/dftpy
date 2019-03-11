from __future__ import print_function
import sys
import numpy as np

from dftpy import structure


def get_struct(s=5.2):
  """Get a structure for silicon"""
  a123 = np.array([[1,1,0],[1,0,1],[0,1,1]])*s # silicon lattice vectors
  atoms = [["Si",[0.,0.,0.]]] # first atom
  atoms += [["Si",[0.25,0.25,0.25]]] # second atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # silicon structure

from dftpy import calculations # calculations module

# quantum espresso
c1 = calculations.Calculation(st,code="QE") # create a DFT calculation
ks1,es1 = c1.bandstructure() # get bandstructure
ks1 = ks1/max(ks1)


# Elk
c2 = calculations.Calculation(st,code="Elk") # create a DFT calculation
ks2,es2 = c2.bandstructure() # get bandstructure
ks2 = ks2/max(ks2)


import matplotlib.pyplot as plt

plt.scatter(ks1,es1,c="red",label="QE")
plt.scatter(ks2,es2,c="blue",label="Elk")

plt.legend()
plt.show()

