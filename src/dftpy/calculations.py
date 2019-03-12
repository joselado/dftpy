from __future__ import print_function
import os
from . import structure
import numpy as np
from . import klist
from . import elkio
from . import qeio
from copy import deepcopy

cores = 5

dftpath = os.path.dirname(os.path.realpath(__file__)) 

elkpath = dftpath+"/codes/elk/src" # path for Elk executable
qepath = dftpath+"/codes/qe/bin" # path for QE executable

h2ev = 13.6*2 # hartree to eV

class Calculation():
  def __init__(self,struct,code="Elk",soc=False,collinear=False,
                   spin=False,m=[0.,0.,0.]): 
    self.code = code # DFT code used
    self.structure = struct # store structure
    self.structure.set_magnetism(m) # add the magnetism
    self.scf = False # scf has not been performed
    self.xc = "PBEsol" # exchange correlation
    self.originalpath = os.getcwd() # go to that path
    self.path = os.getcwd()+"/dfttmp/" # path for the calculation
    self.soc = soc # include spin-orbit coupling
    self.collinear = collinear # collinear calculation
    self.spin = spin # spin polarized calculation
    self.accuracy = 7. # accuracy level from 1 to 10
    self.kmesh = [1,1,1] # kmesh
    self.fermi = 0.0 # Fermi energy
    self.spiral = False # spin spiral calculation
    self.qspiral = [0.,0.,0.] # qvector for spin spiral
    self.energyerror = 1e-5 # error in total SCF energy (in eV)
    self.options = dict() # create dictionary
    os.system("mkdir -p dfttmp") # create the directory
  def get_options(self): return get_options(self)
  def copy(self): return deepcopy(self)
  def dos(self,nk=10):
      if not self.scf: self.gs_energy()
      if self.code=="Elk":
          from .doselk import dos_elk
          return dos_elk(self,nk=nk)
      elif self.code=="QE":
          from .dosqe import dos_qe
          return dos_qe(self,nk=nk)
      else: raise
  def pdos(self,nk=10):
      if not self.scf: self.gs_energy()
      if self.code=="QE":
          from .dosqe import pdos_qe
          return pdos_qe(self,nk=nk)
  def gs_energy(self,restart=False):
    """Get the ground state energy"""
    options = self.get_options() # get additional options
    if self.code=="Elk": # Elk calculation
      if restart:
        structure.write_elk(self.structure,task=1,options=options,
                                path=self.path) 
      else:
        structure.write_elk(self.structure,task=0,options=options,
                                path=self.path) 
      run_elk(self) # run the calculation
      return self.execute(lambda : np.genfromtxt("TOTENERGY.OUT")[-1]*h2ev)
    elif self.code=="QE": # quantum espresso calculation
      structure.write_qe(self.structure,task="scf",options=options,
                                path=self.path) 
      run_qe(self) # run the calculation
      self.fermi = self.execute(lambda: qeio.get_fermi())
      return qeio.get_total_energy(self.path)
    self.scf = True # scf has been performed
  def relax(self,mode="relax"):
    """Perform a structure relaxation"""
    options = self.get_options() # get additional options
    structure.write_qe(self.structure,options=options,
                           path=self.path,task=mode)
    run_qe(self) # run QE
    structure.update_qe_structure(self.structure,self.path)
  def bandstructure(self,**kwargs): 
    return bandstructure(self,**kwargs)
  def execute(self,fun): return execute(self,fun) # execute a function for elk
  def runnscf(self,k=[0.,0.,0.]): runnscf(self,k=k) # non SCF calculation
  def set_spiral(self,q=[0.,0.,0.]):
    self.collinear = False
    self.spiral = True
    self.qspiral = q
    self.spin = True



def run_elk(c):
  """Perform an Elk calculation"""
#  os.system("export OMP_NUM_THREADS="+str(cores)) # cores
  os.chdir(c.path) # go to that path
  os.system(elkpath+"/elk elk.in > elk.info") # run the calculation
  os.chdir(c.originalpath)
#  os.system("export OMP_NUM_THREADS="+str(cores)) # cores



def run_qe(c):
  """Perform an Elk calculation"""
#  os.system("export OMP_NUM_THREADS="+str(cores)) # cores
  os.chdir(c.path) # go to that path
  os.system(qepath+"/pw.x < qe.in > qe.out") # run the calculation
  os.chdir(c.originalpath)
#  os.system("export OMP_NUM_THREADS="+str(cores)) # cores




def execute(c,fun):
  """Execute a command in the Elk folder"""
  os.chdir(c.path) # go to that path
  a = fun()
  os.chdir(c.originalpath)
  return a



def get_options(c):
  """Return a dictionary with different options"""
  opt = dict() # options
  opt["qspiral"] = c.qspiral # spiral calculation
  if c.soc: opt["soc"] = True
  if c.collinear: opt["collinear"] = True
  opt["spin"] = c.spin # spin polarized calculation
  opt["rgkmax"] = c.accuracy
  opt["kmesh"] = c.kmesh
  opt["elkpath"] = elkpath
  opt["spiral"] = c.spiral
  if c.code=="Elk":
      from . import elkinput
      opt["append"] = elkinput.generate_input(c.options) # additional
      opt["append"] += elkinput.extract_attributes(c) # additional
  return opt


def runnscf(self,k=[0.,0.,0.]):
  """Run a non selfconsistent calculation"""
  options = self.get_options() # get additional options
  options["kpoint"] = k # kpoint of the calculation
  options["single_shot"] = True # single shot calculation
  options["single_kpoint"] = True # single shot calculation
  structure.write_elk(self.structure,task=1,options=options,
                              path=self.path) 
  run_elk(self)


def run_lsj_kst(self,states=[1],k=[0.,0.,0.]):
  """Run a non selfconsistent calculation"""
  options = self.get_options() # get additional options
  options["kpoint"] = k # list with the kpoints
  options["single_kpoint"] = True # list with the kpoints
#  options["kpoint"] = k # kpoint of the calculation
#  options["single_shot"] = True # single shot calculation
#  options["single_kpoint"] = True # single shot calculation
  kst = ""
  for s in states: kst += "1   "+str(s)+"\n"
  options["kstlist"] = kst # single shot calculation
  structure.write_elk(self.structure,task=16,options=options,
                              path=self.path) 
  run_elk(self)




def bandstructure(self,**kwargs):
  if not self.scf: self.gs_energy()
  if self.code=="Elk":
    return bandstructure_elk(self,**kwargs)
  elif self.code=="QE":
    return bandstructure_qe(self,**kwargs)
  else: raise



from .bandelk import bandstructure_elk
from .bandqe import bandstructure_qe


