import numpy as np
from . import klist
from . import structure
from . import qeio
from .calculations import run_elk,elkpath,h2ev
import os



def dos_elk(self,nk=None):
    """Return the density of states"""
    options = self.get_options()
    if nk is not None: # if kpoints are given on input
      options["kmesh"] = [nk,nk,nk]
      options["maxscl"] = 1
      options["nempty"] = str(20*len(self.structure.atoms))
      structure.write_elk(self.structure,task=1,options=options,
                                  path=self.path)
      run_elk(self) # diagonalize
    structure.write_elk(self.structure,task=10,options=options,
                                path=self.path)
    run_elk(self) # compute DOS
    m = self.execute(lambda : np.genfromtxt("IDOS.OUT").transpose())
    return m[0]*h2ev,m[1]/h2ev


