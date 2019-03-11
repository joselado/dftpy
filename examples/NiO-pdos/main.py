from __future__ import print_function
import sys
sys.path.append("../../src") # library
import numpy as np

from dftpy import structure



def get_struct(s=7.9):
  """Get a structure for silicon"""
  a123 = np.array([[1,.5,.5],[.5,1.,.5],[.5,.5,1.]])*s # lattice vectors
  atoms = [["Ni",[0.,0.,0.]]] # atom
  atoms += [["Ni",[0.5,0.5,0.5]]] # atom
  atoms += [["O",[0.25,0.25,0.25]]] # atom
  atoms += [["O",[0.75,0.75,0.75]]] # atom
  st = structure.create_structure(a123,atoms) # create structure
  return st



st = get_struct() # NiO structure

from dftpy import calculations # calculations module

# quantum espresso
import matplotlib.pyplot as plt
c1 = calculations.Calculation(st,code="QE",spin=True,m=[0.,0.,1.]) # create a DFT calculation
#es1,ds1 = 
from dftpy import dosqe
pd = c1.pdos(nk=4) # get PDOS object
#exit()
pd = dosqe.read_pdos(c1)
es1,ds1 = pd.energies,pd.pdos_specie[("Ni","up")]
es2,ds2 = pd.energies,pd.pdos_specie[("Ni","dn")]
es1,ds1 = pd.energies,pd.dos_up
es2,ds2 = pd.energies,pd.dos_dn
#es2,ds2 = pd.energies,pd.pdos_specie["O"]
plt.plot(es1,ds1,label="Ni")
plt.plot(es2,ds2,label="O")
plt.legend()
plt.xlim([-5,5])

plt.show()

