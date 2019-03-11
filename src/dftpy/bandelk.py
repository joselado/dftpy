import numpy as np
from . import klist
from . import structure
from . import qeio
from .calculations import run_elk,elkpath,h2ev
import os


def bandstructure_elk(self,kpath=None,spin=False,nk=100):
  """Calculate the bandstructure, together with certain expectation value"""
  self.structure.set_magnetism([0.,0.,0.]) # remove magnetism
  if kpath is None:
    kpath = klist.Kpath().kpoints
  options = self.get_options() # get additional options
  if True: # convetional bandstructure calculation
    options["kinitial"] = kpath[0] # special option
    options["kfinal"] = kpath[-1] # special option
    options["nk"] = nk # special option
    structure.write_elk(self.structure,task=20,options=options,
                            path=self.path)
    run_elk(self) # perform the calculation
    m = self.execute(lambda : np.genfromtxt("BAND.OUT"))
    m[:,1] *= h2ev
    np.savetxt("BANDS.OUT",m)
    return m[:,0],m[:,1]
  else:
    options["bands"] = "1by1" # special option
    fb = open("BANDS.OUT","w") # bandstructure file
    fk = open("KPOINTS.OUT","w") # kpoints file
    for ik in range(len(kpath)): # loop over kpoints
      k = kpath[ik]
      fb.write("\n") # next line
      ik += 1 # increase
      print("Doing ",k)
      options["kpoint"] = k # list with the kpoints
      structure.write_elk(self.structure,task=20,options=options,
                            path=self.path)
      run_elk(self) # perform the calculation
      e = self.execute(lambda : np.genfromtxt("BAND.OUT").transpose()[1])
      for ie in e:
        fb.write(str(ik)+"  ") # write kpoint
        fb.write(str(ie*h2ev)+" \n") # write energy
        fk.write(str(ik)+"  "+str(k[0])+"  "+str(k[1])+"  "+str(k[2])+"\n")
      fb.flush()
      fk.flush()
    fb.close()
    fk.close()
    if spin: # if the spin must be computed
      fs = open("LSJ_BANDS.OUT","w")
      ik = 0 # 
      for ik in range(len(kpath)): # loop over kpoints
        k = kpath[ik] # get the kpoint
        print("Doing ",k)
        self.runnscf(k=k) # run a nonselfconsistent calculation 
        states = range(1,len(e)+1) # as many states
        run_lsj_kst(self,states=states,k=k) # run for states
        for state in states: # loop
          lsj = elkio.read_lsj_kst(self.path,1,state=state)
          fs.write(str(ik)+"  ") # write energy
          fs.write(str(lsj[1][2])+"\n") # write Sz
        fs.flush()
    mout = np.genfromtxt("BANDS.OUT").transpose() # output
    if spin: return mout[0],mout[1],np.genfromtxt("LSJ_BANDS.OUT").transpose()[1]
    else: return mout[0],mout[1]

