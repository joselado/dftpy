from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

import structure


def get_struct(s=5.6):
  """Get a structure for silicon"""
  a123 = np.array([[1,1,0],[1,0,1],[0,1,1]])*s # silicon lattice vectors
  atoms = [["Ga",[0.,0.,0.]]] # first atom
  atoms += [["As",[0.25,0.25,0.25]]] # second atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # silicon structure

import calculations # calculations module
c = calculations.Calculation(st,soc=True) # create a DFT calculation
c.gs_energy() # get ground state energy
kp,es,sz = c.bandstructure(spin=True) # write bandstructure

np.savetxt("BANDS_SZ.OUT",np.matrix([kp,es,sz]).T)
