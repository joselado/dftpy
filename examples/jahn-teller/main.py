from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

import structure


def get_struct(s=8.107):
  """Get a structure for RhTiO3"""
  a123 = np.array([[1,0,0],[0,1,0],[0,0,1]])*s # lattice vectors
  atoms = [["Fe",[0.,0.,0.]]] # atom
  atoms += [["Sr",[0.5,0.5,0.5]]] # atom
  atoms += [["O",[0.5,0.0,0.0]]] # atom
  atoms += [["O",[0.0,0.5,0.0]]] # atom
  atoms += [["O",[0.0,0.0,0.5]]] # atom
  st = structure.create_structure(a123,atoms) # create structure
  return st


import calculations # calculations module

st = get_struct() # get the structure
c = calculations.Calculation(st,m={("Fe",0):[0.,0.,1.]},spin=True) # create a DFT calculation

import optimize

c.code = "QE" # use quantum espresso
c.kmesh = [1,1,1] # gamma point
#c = optimize.volume(c) # optimize the volume
c0 = c.copy() # copy

xs = np.linspace(0.9,1.0,10)
es = []
for x in xs:
  c = c0.copy()
  c.structure.expand([1.,1.,x])
  e = c.gs_energy()
  print(x,e)
  es.append(e)

import matplotlib.pyplot as plt

plt.plot(xs,es,marker="o")
plt.show()

exit()


c = optimize.lattice(c,v=[2]) # optimize the last vector of the lattice
