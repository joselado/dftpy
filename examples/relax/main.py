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


from dftpy import calculations # calculations module

st = get_struct()
c = calculations.Calculation(st,spin=True) # create a DFT calculation
c.relax()

exit()
aa = np.linspace(4.5,5.8,10) # array for the lattice constants
es = [] # array for the energies
for a in aa:
  st = get_struct(a) # get the structure
  c = calculations.Calculation(st) # create a DFT calculation
  e = c.gs_energy() # get the ground state energy
  es.append(e) # store energy
  print("Lattice constant = ",a,"Energy = ",e)

es = np.array(es)
es = es - min(es) # shift energy

import matplotlib.pyplot as plt
plt.plot(aa,es)

plt.show()

