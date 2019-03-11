from __future__ import print_function
import numpy as np
from . import qeio # input and output for Quantum Espresso
from . import elkio # input and output for Elk
from copy import deepcopy # deepcopy
import collections


class Atom():
  """Class for the atom"""
  def __init__(self,name="",r=np.array([0.,0.,0.])):
    self.name = name
    self.r = np.array(r) # position in rfactional coordinates
    self.x = r[0]
    self.y = r[1]
    self.z = r[2]
    self.typeindex = 0 # index of this type
    self.index = 0 # index of this atom
    self.m = [0.,0.,0.]
  def round(self,n=5):
    self.r = np.round(self.r,n)
    self.x = self.r[0]
    self.y = self.r[1]
    self.z = self.r[2]
  def __str__(self):
    out =  "Atom "+self.name+", r = "
    out += str(round(self.r[0],3))+" "
    out += str(round(self.r[1],3))+" "
    out += str(round(self.r[2],3))+" "
    return out
  def __repr__(self): return self.__str__()




class Structure():
  """Structure object"""
  def __init__(self):
    self.dimensionality = 3 # dimensionality
    self.a1 = np.array([1.,0.,0.])
    self.a2 = np.array([0.,1.,0.])
    self.a3 = np.array([0.,0.,1.])
    self.atoms = [] # empty list with the atoms
    self.species = [] # species
  def get_species(self):
    self.species = [] # initialize
    for a in self.atoms:
      if a.name not in self.species: self.species.append(a.name)
    atsp = dict() # dictionary
    for s in self.species: # loop over species
      ii = 0 # initialize
      atsp[s] = [] # empty list
      for a in self.atoms: # loop over atoms
        if a.name==s: 
          a.typeindex = ii
          ii += 1 # increase counter
          atsp[s].append(a) # store position
    self.atomspecies = atsp # store
  def set_magnetism(self,m):
    for a in self.atoms: 
      if type(m)==dict: 
        try: a.m = m[(a.name,a.typeindex)] # set the magnetization
        except: a.m = [0.,0.,0.]
      else: a.m = np.array(m) # set the magnetism
  def round(self,n=4):
    """Round values"""
    self.a1 = np.round(self.a1,n)
    self.a2 = np.round(self.a2,n)
    self.a3 = np.round(self.a3,n)
    for a in self.atoms: a.round(n=n)
    self.get_species()
  def magnetic_inequivalent(self):
    magnetic_inequivalent(self)
  def copy(self):
    return deepcopy(self)
  def expand(self,s=1.0):
    if isinstance(s, collections.Iterable): ss = s # iterable
    else: ss = [s,s,s]  # single number
    self.a1 *= ss[0]
    self.a2 *= ss[1]
    self.a3 *= ss[2]



def read_qe(name="qe.in"):
  """Read a Quantum Espresso structure"""
  st = Structure() # generate a structure object
  aa = qeio.read_lattice(name) # read the lattice vectors
  st.a1,st.a2,st.a3 = aa[0],aa[1],aa[2] # get the lattice vectors
  atoms = qeio.read_atoms(name) # read the atoms
  st.atoms = [Atom(name=a[0],r=a[1]) for a  in atoms]
  st.get_species()
  # get the atoms
  return st # return structure


def update_qe_structure(st,path):
  aa,atoms = qeio.final_structure(path+"qe.out")
  if len(aa)!=0: # lattice vectors optimized
    st.a1,st.a2,st.a3 = aa[0],aa[1],aa[2] # get the lattice vectors
  st.atoms = [Atom(name=a[0],r=a[1]) for a  in atoms]
  st.get_species()
  return st




def write_elk(st,task=0,options=dict(),path=""):
  """Write the Elk input"""
  dictin = dict() # dictionary
  dictin["structure"] = st # structure 
  dictin["task"] = task # task
  from .calculations import elkpath
  dictin["elkpath"] = elkpath # task
  for key in options: dictin[key] = options[key] # Add the options
  elkio.write_elk(dictin,path=path) # write Elk



def write_qe(st0,options=dict(),path="",task="scf"):
  """Write the Elk input"""
  st = st0.copy() # copy
  st.magnetic_inequivalent() # inequivalent
  dictin = dict() # dictionary
  dictin["structure"] = st # structure 
  dictin["task"] = task # structure 
  for key in options: dictin[key] = options[key] # Add the options
  qeio.write_qe(dictin,path=path) # write Elk



def create_structure(a123,atoms):
  """Create a structure"""
  st = Structure() # generate a structure object
  st.a1 = a123[0]
  st.a2 = a123[1]
  st.a3 = a123[2]
  st.atoms = [Atom(name=a[0],r=a[1]) for a in atoms]
  for i in range(len(st.atoms)): st.atoms[i].index = i
  st.get_species() # get the species
  return st



def magnetic_inequivalent(self):
  """If two atoms have non zero initial magnetization,
  assume that they are inequivalent"""
  names = []
  i = 0 # initialize
  for a in self.atoms: # loop over atoms
    m = np.array(a.m)
    if m.dot(m)>0.001: # finite magnetization
      a.name += str(i)
      i += 1 # increase counter
    if a.name not in names: names.append(a.name)
  self.get_species()


