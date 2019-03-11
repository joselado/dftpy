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

c = calculations.Calculation(st,spin=True,m=[0.,0.,-1.]) 
(es,ds) = c.dos()

plt.plot(es,ds)
plt.show()
exit()


c.code = "QE"
c.gs_energy() # get ground state energy
#exit()
c.scf = True
#c.runnscf(k=[0.5,0.5,0.5]) # run a non selfconsistent calculation
kp,es = c.bandstructure(spin=True) # write bandstructure

#sz = sz/np.max(np.abs(sz))
np.savetxt("BAND.OUT",np.matrix([kp,es]).T)
plt.scatter(kp,es)

plt.show()

