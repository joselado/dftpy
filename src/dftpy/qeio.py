
import numpy as np
import os

from .calculations import dftpath
pseudodir = dftpath+"/codes/pseudopotentials"

# dictionary with different cutoff for different elements
cutwf = {}


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
    path = pseudodir
    command = 'ls '+path+'/'+el+".* > .pseudo_files "  # command
    sh(command)  # print in file
    lines = open(".pseudo_files","r").readlines()
    name = el+"." 
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



def write_qe(dictin,path=""):
  """Write the quantum espresso input"""
  def check(name): return (name in dictin)
  names = dictin["structure"].species # list of different atoms
  listat = dictin["structure"].atoms # atoms
  a1 = dictin["structure"].a1 # vector
  a2 = dictin["structure"].a2 # vector
  a3 = dictin["structure"].a3 # vector
  rlat = [a1,a2,a3] # lattice vectors
  fp = '{0:.8f}'.format
  fqe = open(path+'qe.in','w')
  def w(a):
    fqe.write(a+"\n")
  w("&CONTROL")
  w("!                 restart_mode = 'restart' ,")
  w("                 calculation = '"+dictin["task"]+"',")
#  w("                 tstress = .true. ,")
  w("                 verbosity = 'high' ,")
  w("                   wf_collect = .TRUE. ,")
  w("                  pseudo_dir='"+pseudodir+"'")
  w("/")
  w(" &SYSTEM")
  if dictin["task"]=="bands": # bandstructure calculation
    nbnd = sum([pseudo_electrons(get_pseudo(a)) for a in names])
    nbnd = int(nbnd)*3 # convert to integer
    w("     nbnd ="+str(nbnd)+",") # number of bands
  w("                       ibrav = 0")
  scale = 1.0
  w("                   celldm(1) = " + str(scale) + " ,")
  nat = len(listat)
  w("                         nat = "+str(nat)+" ,")
#  names = [a.name for a in listat]
  ntyp = len(names)
  w("                         ntyp = "+str(ntyp)+" ,")
  out = [pseudo_wave_cutoff(get_pseudo(a)) for a in names]
  cutw = np.max([pseudo_wave_cutoff(get_pseudo(a)) for a in names])
  cutd = np.max([pseudo_density_cutoff(get_pseudo(a)) for a in names])
  w("                     ecutwfc = "+str(cutw)+",")
  w("                     ecutrho = "+str(cutd)+",")
  w("                 occupations = 'smearing' ,")
  w("                     degauss = 0.001 ,")
  if check("spin"): # if it has spin
    if dictin["spin"]: # if it is spinful
      w("                    nspin = 2,") # spinful calculation
      for i in range(len(names)):
        m = dictin["structure"].atomspecies[names[i]][0].m[2]
        if abs(m)>1.0: m = m/abs(m)
        w("              starting_magnetization("+str(i+1)+") = "+str(m)+",")
#  w("!                     lda_plus_u = .true. ,")
#  w("!                     lda_plus_u_kind = 0 ,")
#  w("!                     Hubbard_U(1) = 7.0 ,")
  w("/")
  w(" &ELECTRONS")
  w("          startingwfc = 'atomic'")
  # options for good diagonalization
#  w(" !                                diagonalization = 'cg' ,\n!                         diago_full_acc = .TRUE.")
  w("                 mixing_mode = 'TF' ,")
  w("                 mixing_beta = 0.3 ,")
  w("/")
  w(" &IONS")
#  w("!            ion_dynamics = 'damp'\n/")
  w("/")
  w(" &CELL")
  w("    cell_dynamics = 'damp-w'\n")
  w("/")
  w("CELL_PARAMETERS alat")  # print cell parameters
  for i in range(3):
    w("  "+fp(rlat[i][0])+"  "+fp(rlat[i][1])+"  "+fp(rlat[i][2]))
  w("ATOMIC_SPECIES")
  for iname in names:
    w("   "+iname+"  1.000000   "+get_pseudo(iname))
  w("ATOMIC_POSITIONS crystal")
  for iat in listat:
    w("  "+iat.name+"  "+fp(iat.x)+"  "+fp(iat.y)+"  "+fp(iat.z))
  # write the kpoints
  if check("kmesh"): kg = [str(int(ik)) for ik in dictin["kmesh"]]
  else: kg = ["3","3","3"]
  # write the kpath
  if dictin["task"] in ["scf","nscf","relax"]: # selfconsistent calculation
#    if np.max(dictin["kmesh"])==1: w("K_POINTS gamma")
#    else:
      w("K_POINTS automatic")
      w("   "+kg[0]+" "+kg[1]+" "+kg[2]+"   1 1 1")
      w("")
  elif dictin["task"]=="bands": # bandstructure calculation
      w("K_POINTS crystal_b")
      k0 = dictin["kinitial"]
      k1 = dictin["kfinal"]
      nk = dictin["nk"]
      fqe.write("2\n")
      for ik in k0: fqe.write(str(ik)+"  ")
      fqe.write(str(nk)+"\n")
      for ik in k1: fqe.write(str(ik)+"  ")
      fqe.write(str(nk)+"\n")
  fqe.close() # close the file



# get positions of atoms
def read_atoms(name):
  f = open(name).readlines() # get the different atoms
# number of atoms
  nat=int(value(f,'nat'))
# list of atoms
  listat=[atom() for i in range(nat)]
# lines with atoms
  lines=afterpat(f,'ATOMIC_POSITIONS',nat)
  ii=0
  out = [] # output
  for iat in listat:
# get name and coordinates
    line=lines[ii]
    words=line.split()
    at=words[0]
    x=float(words[1])
    y=float(words[2])
    z=float(words[3])
    out.append((at,np.array([x,y,z])))
    ii=ii+1
  return out



# get cell parameters
def read_lattice(name):
  f = open(name).readlines() # read the file
  try: ctyp = int(value(f,"ibrav"))
  except: ctyp = 0
  rlat=np.array([[0.0 for i in range(3)] for j in range(3)])
  if ctyp == 0:    # free cell
    lines=afterpat(f,'CELL_PARAMETERS',3)
    for i in range(3):
      line=lines[i].split()
      for j in range(3):
        rlat[i][j]=float(line[j])
  if ctyp == 8:    # alejandro's cell
    print("Orthorhombic P cell, creating avec")
    a = float(value(f,"celldm(1)"))
    b = float(value(f,"celldm(2)"))
    b = b*a
    c = float(value(f,"celldm(3)"))
    c = c*a
    rlat[0][0] = a
    rlat[1][1] = b
    rlat[2][2] = c
  return rlat



def final_structure(path):
  """Return the final structure"""
  ls = open(path).readlines() # read the lines
  name = "/tmp/qe.final"
  found = False
  f = open(name,"w")
  for l in ls: # loop over lines
    if "Begin final coordinates" in l:
      found = True
    if found: f.write(l)
    if "End final coordinates" in l:
      break
  f.close()
  # get the lattice vector
  ls = open(name).readlines() 
  found = False 
  c = 0 # counter
  aa = []
  for l in ls:
    if found:
      c += 1 # increase counter
      a = l.split() # split line
      aa.append(np.array([float(a[i]) for i in range(3)]))
    if c==3: break # break
    if "CELL_PARAMETERS" in l:
      found = True
  # read the atoms
  ls = open(name).readlines() 
  out = [] # list of atoms
  found = False
  for l in ls:
    if "End final coordinates" in l: break
    if found:
      a = l.split() # split line
      out.append([a[0],[float(a[1]),float(a[2]),float(a[3])]])
    if "ATOMIC_POSITIONS" in l:
      found = True
  return aa,out



def get_total_energy(path):
  """Get the total energy"""
  ls = open(path+"qe.out").readlines()
  for l in ls:
    if "!    total energy" in l:
      t = l.split("=")[1].split()[0]
      e = float(t)*13.6 # to eV
  return e



def pseudo_wave_cutoff(pseudo):
  """Get the energy cutoff of a certain pseudopotential"""
  ls = open(pseudodir+"/"+pseudo).readlines() 
  i = 0
  for l in ls:
    i += 1
    if "Suggested minimum cutoff for wavefunctions" in l:
      out = l.split(":")[1].split("Ry")[0]
      return float(out)
    if i>30: break
  return 50.0



def pseudo_electrons(pseudo):
  """Get the number of electrons a certain pseudopotential"""
  ls = open(pseudodir+"/"+pseudo).readlines() 
  i = 0
  for l in ls:
    i += 1
    if "z_valence" in l:
      out = l.split("=")[1].replace('"','')
      return float(out)
    if i>200: break
  print("Electrons not found in",pseudo)
  return 10.0


def pseudo_density_cutoff(pseudo):
  """Get the energy cutoff of a certain pseudopotential"""





def pseudo_density_cutoff(pseudo):
  """Get the energy cutoff of a certain pseudopotential"""
  ls = open(pseudodir+"/"+pseudo).readlines() 
  i = 0
  for l in ls:
    i += 1
    if "Suggested minimum cutoff for charge density" in l:
      out = l.split(":")[1].split("Ry")[0]
      return float(out)
    if i>30: break
  return 50.0





def write_bandsin(spin=1):
    """Write bands.in"""
    f = open("bands.in","w")
    f.write("&BANDS\n")
    f.write("filband = 'band_QE.dat' ,\n")
    f.write("spin_component = "+str(spin)+" ,\n /\n")
    f.close()



def write_dosin(delta=0.5):
    """Write bands.in"""
    f = open("dos.in","w")
    f.write("&DOS\n")
    f.write("fildos = 'dos.dat' ,\n")
    f.write("degauss = "+str(delta/13.6)+" ,\n")
    f.write("deltaE = "+str(delta/13.6/4.0)+" ,\n")
    f.write("/\n")
    f.close()


def write_pdosin(delta=0.1):
    """Write bands.in"""
    f = open("projwfc.in","w")
    f.write("&projwfc\n")
#    f.write("fildos = 'dos.dat' ,\n")
    f.write("degauss = "+str(delta/13.6)+" ,\n")
#    f.write("deltaE = "+str(delta/13.6/4.0)+" ,\n")
    f.write("/\n")
    f.close()









def get_fermi():
    name = "qe.out"
    ls = open(name).readlines()
    fermi= None
    for l in ls:
        if "the Fermi energy is" in l:
            fermi = float(l.split("is")[1].split("ev")[0])
    if fermi is None: 
        print("Warning, fermi energy not found in",name)
    return fermi
