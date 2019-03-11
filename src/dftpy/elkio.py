
import os
import numpy as np

def get_element(name):
  """Gets the true name of an element, assumes that there might be a number"""
  out = "" # initialize
  if name[-1] in [str(i) for i in range(10)]: 
    for i in range(len(name)-1): # except the last one
      out += name[i] # add to the string
    return get_element(out) # recall the function, just in case two digits
  return name


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
def read_atoms(name):
  f = read(name).readlines()
  lines=afterpat(f,'atoms',-1)
  nty = int(lines[0].split()[0]) # number of types
  ili = 1
  ity = 1
  out = []
  for ity in range(nty):
    name = lines[ili].split("'")[1].split(".")[0] # read name
    ili += 1   # next line
    inat = int(lines[ili].split()[0]) # read number of this type
    ili += 1   # next line  
    for ia in range(inat):
      words=lines[ili].split()
      ili += 1
      iat = atom()
      name = name
      x = float(words[0])
      y = float(words[1])
      z = float(words[2])
      out.append((name,np.array([x,y,z])))
  return out

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
  return (dnames,ntyp)





def write_elk(dictin,path=""):
  """Write Elk input"""
  def check(name): return (name in dictin)
  task = dictin["task"] # calculation to perform
  nameats = dictin["structure"].species # different species
  listat = dictin["structure"].atoms # atoms
  a1 = dictin["structure"].a1 # vector
  a2 = dictin["structure"].a2 # vector
  a3 = dictin["structure"].a3 # vector
  rlat = [a1,a2,a3] # lattice vectors
  ntyps = [len(dictin["structure"].atomspecies[s]) for s in nameats]

  fp = '{0:.8f}'.format

  felk = open(path+'elk.in','w') # elk input
  def w(a,endl=True):
    felk.write(a)
    if endl: felk.write("\n")

  def checkwrite(name):
      """Write a variable if it exist in the dictionary"""
      if check(name):  w(name+"\n"+obj2elk(dictin[name])+"\n\n") 

  w("tasks\n"+str(task)+"\n\n")
  w("")
  w("")
  w("! GGA exchange correlation")
  w("xctype")
  w("  20  0  0\n")  # LDA potential
  if check("collinear"): # collinear calculation
    if dictin["collinear"]: w("cmagz\n.true.\n") # use collinear formalism
  # spin polarized calculation
  w("avec")
  w("    "+fp(rlat[0][0])+"  "+fp(rlat[0][1])+"   "+fp(rlat[0][2]))
  w("    "+fp(rlat[1][0])+"  "+fp(rlat[1][1])+"   "+fp(rlat[1][2]))
  w("    "+fp(rlat[2][0])+"  "+fp(rlat[2][1])+"   "+fp(rlat[2][2]))
  w("")
#  path = get_elkpath() # get the elk path
  path = dictin["elkpath"] # get the path
  checkwrite("nempty") # empty states
  w("sppath\n  '"+path+"/../species/'\n")
  w("atoms")
  w("  "+str(len(nameats))+"              : nspecies")
  for i in range(len(nameats)):
    w("  "+"'"+nameats[i]+".in'                  : spfname")
    w("  "+str(ntyps[i])+"                      : natoms")
    for iat in listat:
      if iat.name == nameats[i]:
        w("  "+fp(iat.x)+"  "+fp(iat.y)+"  "+fp(iat.z)+"   ",endl=False)
        w("  "+fp(iat.m[0])+"  "+fp(iat.m[1])+"  "+fp(iat.m[2])+"   ")
  w("")
  w("")
  # write the kmesh
  if check("kmesh"): kp = [str(int(ik)) for ik in dictin["kmesh"]]
  else: kp = ["3","3","3"]
  w("ngridk")
  w("  "+str(kp[0])+" "+str(kp[1])+" "+str(kp[2]))
  w("")
  w("")
  if task in [20]: # 1d band structure plot
    w("plot1d")
    if check("bands"): # do a single k-point
      if dictin["bands"]=="1by1": # do a single k-point
        w("1 1")
        k = dictin["kpoint"] # kpoint
        for ik in k: w(str(ik)+" ",endl=False)
        w("\n\n")
    else: # conventional bandstructure
      w("2  200")
      for ik in dictin["kinitial"]: w(str(ik)+" ",endl=False)
      w("")
      for ik in dictin["kfinal"]: w(str(ik)+" ",endl=False)
      w("\n\n")
    
  # additional parameters
  if check("soc"): w("spinorb\n"+obj2elk(dictin["soc"])+"\n") # spin orbit 
  if check("single_shot"): # if the parameter exists
    if dictin["single_shot"]: w("maxscl\n1\n")
  # check if it is given on input
  if check("maxscl"): w("maxscl\n"+str(dictin["maxscl"])+"\n\n")
  if check("single_kpoint"): # if the parameter exists
    if dictin["single_kpoint"]: # compute a single kpoint
      w("ngridk\n1  1  1\n") # set the mesh
      k = np.array(dictin["kpoint"]) # get the kpoint
      w("vkloff\n"+obj2elk(k)+"\n") # write the shift
  if check("kstlist"): # if the parameter exists
      w("kstlist\n"+obj2elk(dictin["kstlist"])+"\n") # add this parameter
  if check("rgkmax"):
      w("rgkmax\n"+str(dictin["rgkmax"])+"\n")
  else: w("\nrgkmax\n  7.0\n") # RGKMAX
  if check("spin"): # if it has spin degree of freedom
    if dictin["spin"]: 
      w("spinpol\n.true.\n") # spin polarized calculation
      w("reducebf\n0.7\n") # remove exchange fields
  
  # spin spiral calculation
  if check("spiral"):
      w("spinsprl\n.true.\n\n") # spin spiral calculation
      q = dictin["qspiral"] 
      w("vqlss\n"+str(q[0])+"  "+str(q[1])+"  "+str(q[2])+"\n")





def obj2elk(a):
  """Transform an object to the elk format"""
  if type(a)==type(True):
    if a: return ".true."
    else: return ".false."
  elif type(a)==type("a"): return a
  elif type(a)==type(1.2): return str(a)
  elif type(a)==type(np.array([0.,0.,0.])):
    o = ""
    for ia in a: o+=str(ia)+" "
    return o
  else: raise





def get_elkpath():
  """Return the path for Elk"""
  os.system("which elk > /tmp/elkpath.txt")
  path = open("/tmp/elkpath.txt").read()
  path = path.replace("/src/elk\n","") 
  return path 



def read_lsj_kst(path,k=1,state=1,specie=1,atom=1):
  """Read the Elk file LSJ_KST.OUT"""
  l = open(path+"/LSJ_KST.OUT").read() # open the file
  l = l.split("\n") # split in different lines
  ii = 0
  for ii in range(len(l)): # loop over lines
     il = l[ii] # get this line
     if "k-point :" in il: # if kpoint line
       ik = int(float(il.split()[2])) # read the kpoint
       if ik==k: # found the kpoint
         il = l[ii+1] # get the next line
         if int(float(il.split()[2]))==state: # found the state
           il = l[ii+2] # get the next line
           if int(float(il.split()[2]))==specie: # found the specie
             if int(float(il.split()[6]))==atom: # found the atom
               il = l[ii+3] # get the next line
               il = il.split()
               lout = [float(il[2]),float(il[3]),float(il[4])]
               il = l[ii+4] # get the next line
               il = il.split()
               sout = [float(il[2]),float(il[3]),float(il[4])]
               il = l[ii+5] # get the next line
               il = il.split()
               jout = [float(il[2]),float(il[3]),float(il[4])]
               return np.array(lout),np.array(sout),np.array(jout) # return
  raise # error










