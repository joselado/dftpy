from __future__ import print_function
import sys
import numpy as np

from dftpy import structure
from dftpy import calculations # calculations module
import matplotlib.pyplot as plt


def get_struct(s=2.7):
  """Get a structure for Fe"""
  a123 = np.array([[1,1,-1],[1,-1,1],[-1,1,1]])*s # lattice vectors
  atoms = [["Co",[0.,0.,0.]]] # first atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # Fe structure

c = calculations.Calculation(st,spin=True,m=[0.,0.,-1.],code="Elk") 
c.soc = True
c.options["autolinengy"] = True # set linenergy to True
c.options["socscf"] = 4 # set linenergy to True
c.gs_energy() # compute ground state energy



