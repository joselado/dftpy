import numpy as np

from .elk2qe import write_qe

def getats_vasp(ll):
  """Get atoms from the vasp input"""
  class Atom(): pass
  read = False
  names = get_names(ll) # get names of the atoms
  rs = []
  ii = 0 # counter for atoms
  for l in ll:
    if read and ii<len(names):
      l = l.split()
      a = Atom()
      a.name = names[ii]  # append name
      a.x = float(l[0]) 
      a.y = float(l[1]) 
      a.z = float(l[2]) 
      rs.append(a) # append the atom
      ii += 1
    if "Direct" in l:
      read = True
  return rs # return list of atoms


def get_names(ll):
  names = ll[5].split()
  rep = ll[6].split()
  rep = [int(float(i)) for i in rep] # number of repetitions
  no = []
  for (n,r) in zip(names,rep):
    for i in range(r):
      no.append(n) # append name
  return no

def get_rlat_vasp(ll):
  """Read the lattice in the vasp input"""
  scale = 1./0.529  # from bhor to ams
  l = ll[2].split()
  r1 = np.array([float(l[0]),float(l[1]),float(l[2])])*scale
  l = ll[3].split()
  r2 = np.array([float(l[0]),float(l[1]),float(l[2])])*scale
  l = ll[4].split()
  r3 = np.array([float(l[0]),float(l[1]),float(l[2])])*scale
  return [r1,r2,r3]


# hack to use elk2qe
if __name__=="__main__":
  ll = open("vasp.in").readlines()  # read vasp input
  listat = getats_vasp(ll) # get atoms
  rlat = get_rlat_vasp(ll) # get the reciprocal lattice
  write_qe(ll,listat=listat,rlat=rlat)  # write QE input



