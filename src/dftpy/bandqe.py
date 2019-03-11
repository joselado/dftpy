import numpy as np
from . import klist
from . import structure
from . import qeio
from .calculations import run_qe,qepath
import os


def bandstructure_qe(self,kpath=None,spin=False):
  """Calculate the bandstructure, together with certain expectation value"""
  self.structure.set_magnetism([0.,0.,0.]) # remove magnetism
  if kpath is None:
    kpath = klist.Kpath().kpoints
  options = self.get_options() # get additional options
  options["kinitial"] = kpath[0] # special option
  options["kfinal"] = kpath[-1] # special option
  options["nk"] = 100 #len(kpath) # special option
  structure.write_qe(self.structure,task="bands",options=options,
                          path=self.path)
  run_qe(self) # run QE
  if self.spin: # spinful calculation
    self.execute(lambda : qeio.write_bandsin(1)) # write bands.in
    self.execute(lambda : os.system(qepath+"/bands.x < bands.in > bands.out")) # run
    m1 = self.execute(lambda : np.genfromtxt("band_QE.dat.gnu"))
    self.execute(lambda : qeio.write_bandsin(2)) # write bands.in
    self.execute(lambda : os.system(qepath+"/bands.x < bands.in > bands.out")) # run
    m2 = self.execute(lambda : np.genfromtxt("band_QE.dat.gnu"))
    ks = np.concatenate([m1[:,0],m2[:,0]]) # kpoints
    es = np.concatenate([m1[:,1],m2[:,1]]) # energies
    m = np.array([ks,es]).transpose()
  else: # spinless calculation
    self.execute(lambda : qeio.write_bandsin()) # write bands.in
    self.execute(lambda : os.system(qepath+"/bands.x < bands.in > bands.out")) # run
    m = self.execute(lambda : np.genfromtxt("band_QE.dat.gnu"))
  m[:,1] -= self.fermi # minus the fermi energy
  np.savetxt("BANDS.OUT",m)
  return m[:,0],m[:,1]
#  self.execute(lambda : qeio.write_bandsin(self)) # write bands.in
#  run_elk(self) # perform the calculation
#  m = self.execute(lambda : np.genfromtxt("BAND.OUT"))
#  m[:,1] *= h2ev
#  np.savetxt("BANDS.OUT",m)
#  return m[:,0],m[:,1]


