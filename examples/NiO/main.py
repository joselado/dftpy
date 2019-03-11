from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

import structure


def get_struct(s=7.9):
  """Get a structure for silicon"""
  a123 = np.array([[1,.5,.5],[.5,1.,.5],[.5,.5,1.]])*s # lattice vectors
  atoms = [["Ni",[0.,0.,0.]]] # atom
  atoms += [["Ni",[0.5,0.5,0.5]]] # atom
  atoms += [["O",[0.25,0.25,0.25]]] # atom
  atoms += [["O",[0.75,0.75,0.75]]] # atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # silicon structure

import calculations # calculations module
maf = {("Ni",0):[0.,0.,1.],("Ni",1):[0.,0.,-1.]}
c = calculations.Calculation(st,spin=True,m=maf) 
#c.code = "QE"
eaf = c.gs_energy() # get ground state energy
c = calculations.Calculation(st,spin=True,m=[0.,0.,1.]) 
#c.code = "QE"
efe = c.gs_energy() # get ground state energy


print("AF-FE",eaf-efe)
