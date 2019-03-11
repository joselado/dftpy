from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

import structure


def get_struct(s=5.2):
  """Get a structure for silicon"""
  a123 = np.array([[1,1,0],[1,0,1],[0,1,1]])*s # silicon lattice vectors
  atoms = [["Si",[0.,0.,0.]]] # first atom
  atoms += [["Si",[0.25,0.25,0.25]]] # second atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # silicon structure

import calculations # calculations module
c = calculations.Calculation(st) # create a DFT calculation
#c.spin = True
c.code = "QE"
c.gs_energy() # get ground state energy
c.scf = True
ks,es = c.bandstructure() # write bandstructure

kp = ks
#kp = np.array([k[0] for k in ks])
kp = kp/max(kp)

print(kp.shape,es.shape)

import matplotlib.pyplot as plt

plt.scatter(kp,es,c="black")

plt.show()

