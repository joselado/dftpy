import numpy as np
from . import klist
from . import structure
from . import qeio
from .calculations import run_qe,qepath
import os
from .filetk import get_files



def dos_qe(self,nk=None):
    """Return the density of states"""
    options = self.get_options()
    if nk is not None: # if kpoints are given on input
      options["kmesh"] = [nk,nk,nk]
      structure.write_qe(self.structure,task="nscf",options=options,
                                  path=self.path)
      run_qe(self) # diagonalize
    m = self.execute(lambda : qeio.write_dosin()) # write DOS input
    self.execute(lambda : os.system(qepath+"/dos.x < dos.in > dos.out"))
#    structure.write_elk(self.structure,task=10,options=options,
#                                path=self.path)
#    run_elk(self) # compute DOS
    m = self.execute(lambda : np.genfromtxt("dos.dat").transpose())
    return m[0]-self.fermi,m[1]/2. # not sure why this 2 should be here



def pdos_qe(self,nk):
    """Compute the partial DOS, and return an PDOS object"""
    dos_qe(self,nk=nk) # run the DOS
    self.execute(lambda :qeio.write_pdosin()) # write input
    self.execute(lambda : os.system("rm -f pwscf.pdos_*")) # remove old files
    self.execute(lambda : os.system(qepath+"/projwfc.x < projwfc.in > projwfc.out"))
    return read_pdos(self)


def read_pdos(self):
    """Return the PDOS object"""
    from .dos import PartialDos
    pd = PartialDos(self) # PDOS object
    # read the different DOS
    fs = get_files(self.path,".pdos_atm") # get PDOS files
    fermi = self.execute(lambda :qeio.get_fermi()) # get fermi energy
    for f in fs: # loop over PDOS files
#        print(f) # this is just a check
        atn = int(f.split("atm#")[1].split("(")[0]) # atom number
        atom = f.split("(")[1].split(")")[0] # atom number
        l = f.split("(")[2].split(")")[0] # l component
        wf = int(f.split("#")[2].split("(")[0]) # wavefunction number
        m = self.execute(lambda : np.genfromtxt(f).transpose())
        pd.set_energies(m[0]-fermi) # set the energies
        if self.spin: # spinful calculation
          pd.add_pdos(atn-1,l,wf,m[1],spin="up") # set the DOS
          pd.add_pdos(atn-1,l,wf,m[2],spin="dn") # set the DOS
        else: # spinless calculation
          pd.add_pdos(atn-1,l,wf,m[1],spin=None) # set the DOS
    pd.add_contributions() # do some post processing
    m = self.execute(lambda : np.genfromtxt("pwscf.pdos_tot").transpose())
    if self.spin: # spinful calculation
      pd.add_total(m[1],spin="up")
      pd.add_total(m[2],spin="dn")
      pd.add_total(m[1]+m[2],spin="dn")
    else:
      pd.add_total(m[1],spin=None)
    return pd # Return PartialDOS object


