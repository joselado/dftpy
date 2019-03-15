from __future__ import print_function
import sys
import numpy as np

from dftpy import structure


def get_struct(s=4.5):
  """Get a structure for silicon"""
  a123 = np.array([[1,1,0],[1,0,1],[0,1,1]])*s # silicon lattice vectors
  atoms = [["Si",[0.,0.,0.]]] # first atom
  atoms += [["Si",[0.25,0.25,0.25]]] # second atom
  st = structure.create_structure(a123,atoms) # create structure
  return st


from dftpy import calculations # calculations module
from dftpy import optimize

st = get_struct() # get the structure
calculations.cores = 4
c = calculations.Calculation(st,code="Elk") # create a DFT calculation
c = optimize.volume(c) # optimize the volume
