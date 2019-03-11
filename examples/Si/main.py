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


aa = np.linspace(4.5,5.8,10) # array for the lattice constants
es1 = [] # array for the energies
es2 = [] # array for the energies
for a in aa:
  st = get_struct(a) # get the structure
  c = calculations.Calculation(st) # create a DFT calculation
  c.code="QE" # Use quantum Espresso
  e1 = c.gs_energy() # get the ground state energy
  c.code="Elk" # Use Elk
  e2 = c.gs_energy() # get the ground state energy
  es1.append(e1) # store energy
  es2.append(e2) # store energy
  print("Lattice constant = ",a,"Energy = ",e1,e2)

es1 = np.array(es1)
es2 = np.array(es2)
es1 = es1 - min(es1) # shift energy
es2 = es2 - min(es2) # shift energy

import matplotlib.pyplot as plt
plt.plot(aa,es1,marker="o",label="QE")
plt.plot(aa,es2,marker="o",label="Elk")

plt.legend()
plt.show()

