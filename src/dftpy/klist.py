from __future__ import print_function
import numpy as np


class Kpath(): # class for the klist
  def __init__(self):
    self.kpoints = [np.array([1.,1.,1.])*k for k in np.linspace(0,1.,100)]
    labels = ["" for k in self.kpoints] # nothing
    labels[0] = "$\Gamma$"
    labels[-1] = "$\Gamma$"
    self.labels = labels # store








