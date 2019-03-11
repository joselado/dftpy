from dftpy import structure
from dftpy import calculations # calculations module
from dftpy import dosqe # calculations module

st = structure.read_qe("qe.in") # read Quantum Espresso structure

c = calculations.Calculation(st,code="QE") 
#pd = c.pdos(nk=10) # get PDOS object
pd = dosqe.read_pdos(c)
#es,ds = pd.energies,pd.pdos_specie[("Re")]
e0,d0 = pd.energies,pd.get_sum([None,"s",None]) # get the s PDOS
e1,d1 = pd.energies,pd.get_sum([None,"d",None]) # get the d PDOS
e2,d2 = pd.energies,pd.get_sum([None,"p",None]) # get the d PDOS
import matplotlib.pyplot as plt
import numpy as np

np.savetxt("PDOS_S_QE.OUT",np.array([e0,d0]).T)
np.savetxt("PDOS_D_QE.OUT",np.array([e1,d1]).T)
np.savetxt("PDOS_P_QE.OUT",np.array([e2,d2]).T)


plt.plot(e0,d0)
plt.plot(e1,d1)
plt.plot(e2,d2)
plt.show()

