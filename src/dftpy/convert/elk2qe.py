#!/usr/bin/python

import numpy as np
import os


# dictionary with different cutoff for different elements
cutwf = {}
#cutwf{"C"} = (40,400)
#cutwf{"O"} = (40,400)
#cutwf{"H"} = (40,400)


def get_element(name):
  """Gets the true name of an element, assumes that there might be a number"""
  out = "" # initialize
  if name[-1] in [str(i) for i in range(10)]: 
    for i in range(len(name)-1): # except the last one
      out += name[i] # add to the string
    return get_element(out) # recall the function, just in case two digits
  return name

def get_pseudo(el,xc="pbe",rel=False,paw=True):
  """Gets the pseudopotential of a particular element"""
  try:
    el = get_element(el)
    sh = os.system
    cd = os.chdir
    path = '/share/inl/pseudopotentials' 
    command = 'ls '+path+'/'+el+".* > .pseudo_files "  # command
    sh(command)  # print in file
    lines = open(".pseudo_files","r").readlines()
    if not rel: # if not relativistic
      name = el+"."+xc+"-"  
    if rel:  # if relativistic
      name = el+".rel-"+xc+"-"  
    for l in lines:
      l = l.split("/")
      l = l[len(l)-1]
      l = l.split("\n")
      l = l[0]
      if name in l:
        if paw:
          if "paw" in l:
            break
        if not paw:
          if not "paw" in l:
            break
    return l
  except: return ""






# return lines with a pattern
def pattern(f,pat):
  out=[]
  for line in f:
    if pat in line:
      out=out+[line]
  return out

# gets the value of certain parameter
def value(f,var):
  out=pattern(f,var)[0]
  out=out.split(',')
  if len(out)<3: # one parameter per line
    out = out[0]
    out=out.split('=')[1]
    return out
  if len(out)>=3: # several parameters per line
    for ip in out:
      if var in ip:
        ip=ip.split('=')[1]
        return ip
        

# return lines after a pattern
def afterpat(f,pat,nl):
  """ Get nl lines after the pattern pat in lines f"""
  store=False
  listline=[]
  count=0
  for line in f:
    if store:
      listline=listline+[line]
      count=count+1
    if pat in line:
      store=True
    if count == nl:
      break 
  return listline




# atom class
class atom():
  name='Atom'
  x=0.0
  y=0.0
  z=0.0



# get positions of atoms
def getats_elk(f):
  lines=afterpat(f,'atoms',-1)
  nty = int(lines[0].split()[0]) # number of types
  print "Found "+str(nty)+" different species"
  ili = 1
  ity = 1
  listat = []
  for ity in range(nty):
    name = lines[ili].split("'")[1].split(".")[0] # read name
    ili += 1   # next line
    inat = int(lines[ili].split()[0]) # read number of this type
    ili += 1   # next line  
    for ia in range(inat):
      words=lines[ili].split()
      ili += 1
      iat = atom()
      iat.name = name
      iat.x = float(words[0])
      iat.y = float(words[1])
      iat.z = float(words[2])
      listat += [iat]
  return listat

# contract and scale positions
def mvsc(at,mv=[0.0,0.0,0.0],sc=[1.0,1.0,1.0]):
  atout=atom()
  atout.name=at.name
  atout.x=sc[0]*at.x+mv[0]
  atout.y=sc[1]*at.y+mv[1]
  atout.z=sc[2]*at.z+mv[2]
  return atout 


# get cell parameters
def get_rlat_elk(f):
  lines = afterpat(f,"avec",3)
  rlat = [[],[],[]]
  for i in range(3):
    l = lines[i].split()
    rlat[i] = [float(l[0]),float(l[1]),float(l[2])]
  return rlat
 

# get different atom names
def get_atnames(listat):
  dnames = []            # list with name of atoms
  for at in listat:
    if not at.name in dnames:
      dnames = dnames + [at.name]
  ntyp = [0 for i in range(len(dnames))]    # list with number of atoms
  for i in range(len(dnames)):
    for at in listat:
      if at.name == dnames[i]:
        ntyp[i] += 1
    print "Found  "+str(ntyp[i])+ "  "+ dnames[i] 
  return (dnames,ntyp)



###################################
##################################
# function to create qe.in
#################################
#################################




def write_qe(f,listat=None,rlat=None):
  """Write the quantum espresso input"""
  if listat==None:
    listat=getats_elk(f)
  if rlat==None:
    rlat = get_rlat_elk(f)

  fp = '{0:.8f}'.format

  fqe = open('qe.in','w')
  out = "" # empty string
  def w(a): out += a # append
  w("&CONTROL")
  w("!                 restart_mode = 'restart' ,")
  w("                 calculation = 'scf' ,")
  w("                 tstress = .true. ,")
  w("                 verbosity = 'high' ,")
  w("                   wf_collect = .TRUE. ,")
  w("                  pseudo_dir='/share/inl/pseudopotentials/'")
  w("/")
  w(" &SYSTEM")
  w("                       ibrav = 0")
  scale = 1.0
  if "scale" in f:   # if scale in elk.in get it
    scale = float(afterpat(f,"scale",1)[0])
    print "found scale in elk.in",scale
  w("                   celldm(1) = " + str(scale) + " ,")
  nat = len(listat)
  w("                         nat = "+str(nat)+" ,")
  names = get_atnames(listat)[0]
  ntyp = len(names)
  w("                         ntyp = "+str(ntyp)+" ,")
  w("                     ecutwfc = 60 ,")
  w("                     ecutrho = 600 ,")
  w("                 occupations = 'smearing' ,")
  w("                     degauss = 0.002 ,")
  w("!                     nspin = 2 ,")
  w("!                     lda_plus_u = .true. ,")
  w("!                     lda_plus_u_kind = 0 ,")
  w("!                     Hubbard_U(1) = 7.0 ,")
  w("!                     starting_magnetization(1) = 1 ,")
  w("/")
  w(" &ELECTRONS")
  w("          startingwfc = 'atomic'")
  # options for good diagonalization
  w(" !                                diagonalization = 'cg' ,\n!                         diago_full_acc = .TRUE.")
  w("                 mixing_mode = 'local-TF' ,")
  w("/")
  w(" &IONS")
  w("            ion_dynamics = 'damp'\n/")
  w(" &CELL")
  w("    cell_dynamics = 'damp-w'\n")
  w("/")
  w("CELL_PARAMETERS alat")  # print cell parameters
  for i in range(3):
    w("  "+fp(rlat[i][0])+"  "+fp(rlat[i][1])+"  "+fp(rlat[i][2]))
  w("ATOMIC_SPECIES")
  for iname in names:
    w("   "+iname+"  1.000000   "+get_pseudo(iname))
#    w("   "+iname+"  1.000000   "+iname +".pbe-mt_fhi.UPF")
  w("ATOMIC_POSITIONS crystal")
  for iat in listat:
    w("  "+iat.name+"  "+fp(iat.x)+"  "+fp(iat.y)+"  "+fp(iat.z))
  w("K_POINTS automatic")
  try:  kg = afterpat(f,"ngridk",1)[0].split()  # kgrid
  except:  kg = ["3","3","3"]
  w("   "+kg[0]+" "+kg[1]+" "+kg[2]+"   1 1 1")
  w("")
  return out # return string

