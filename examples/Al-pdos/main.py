from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

from dftpy import structure


def get_struct(s=3.83):
  """Get a structure for silicon"""
  a123 = np.array([[1,1,0],[1,0,1],[0,1,1]])*s # silicon lattice vectors
  atoms = [["Si",[0.,0.,0.]]] # first atom
  st = structure.create_structure(a123,atoms) # create structure
  return st

st = get_struct() # silicon structure

from dftpy import calculations # calculations module

# quantum espresso
import matplotlib.pyplot as plt
c1 = calculations.Calculation(st,code="QE") # create a DFT calculation
#es1,ds1 = 
pd = c1.pdos(nk=4) # get bandstructure
#c2 = calculations.Calculation(st,code="Elk") # create a DFT calculation
#es2,ds2 = c2.dos(nk=8) # get bandstructure
es1,ds1 = pd.energies,pd.pdos_atom[0]
plt.plot(es1,ds1,label="QE")
#plt.plot(es2,ds2,label="Elk")
#plt.legend()

plt.show()

