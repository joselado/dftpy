from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

from dftpy import structure


def get_struct(s=2.7):
  """Get a structure for Fe"""
  a123 = np.array([[1,1,-1],[1,-1,1],[-1,1,1]])*s # lattice vectors
  atoms = [["Fe",[0.,0.,0.]]] # first atom
  st = structure.create_structure(a123,atoms) # create structure
  return st


st = get_struct() # Fe structure

from dftpy import calculations # calculations module

# quantum espresso
import matplotlib.pyplot as plt
c1 = calculations.Calculation(st,code="QE",spin=True,m=[0.,0.,1.]) # create a DFT calculation
#c1.gs_energy()
#es1,ds1 = 
from dftpy import dosqe
#pd = c1.pdos(nk=4) # get PDOS object
#exit()
pd = dosqe.read_pdos(c1)
#es1,ds1 = pd.energies,pd.pdos_specie[("Fe","up")]
#es2,ds2 = pd.energies,pd.pdos_specie[("Fe","dn")]
es1,ds1 = pd.energies,pd.dos_up
es2,ds2 = pd.energies,pd.dos_dn
#es2,ds2 = pd.energies,pd.pdos_specie["O"]
plt.plot(es1,-ds1,label="Up")
plt.plot(es2,ds2,label="Dn")
plt.plot(es2,-ds1+ds2,label="Difference",c="black")
plt.legend()
plt.xlim([-5,5])

plt.show()

