import numpy as np

class PartialDos():
    def __init__(self,c):
        self.atoms = c.structure.atoms
        self.structure = c.structure # store the structure
        self.pdos = dict() # dictionary
        self.spin = c.spin # spin degree of freedom
    def set_energies(self,x):
        self.energies = x # initialize energies
    def add_pdos(self,atom,l,wf,y,spin=None):
        """Add a new PDOS to the object"""
        n = self.structure.atomspecies
        if spin is None: self.pdos[(atom,l,wf)] = y.copy() # store this PDOS
        else: self.pdos[(atom,l,wf,spin)] = y.copy() # store this PDOS
    def get_sum(self,ls):
        """
        Return the sum of PDOS according to a list
        """
        out = np.zeros(self.energies.shape[0]) # initialize
        for key in self.pdos: # loop over all the pdos
            add = True
            for (il,ik) in zip(ls,key): # loop over entries
                if il is None: continue # next iteration, accept it
                if str(il)==str(ik): continue # next iteration, accept it
                add = False ; break # this contribution will not be added
            if add: out += self.pdos[key] # add this contribution
        return out # return contribution
    def add_contributions(self): # add the different contributions
        self.pdos_atom = dict() # dictionary
        # contributions for a specific atom
        for i in range(len(self.atoms)): # loop over atoms
            self.pdos_atom[i] = np.zeros(len(self.energies)) # initialize
            if self.spin: # spinful calculation
              self.pdos_atom[(i,"up")] = np.zeros(len(self.energies)) 
              self.pdos_atom[(i,"dn")] = np.zeros(len(self.energies)) 
            for key in self.pdos: # loop over dictionary
                if key[0]==i: # if for that atom
                    self.pdos_atom[i] += self.pdos[key] # add contribution
                    if self.spin: # spinful calculation
                      self.pdos_atom[(i,key[3])] += self.pdos[key] 
        # contributions for a specific specie
        self.pdos_specie = dict() # dictionary
        for s in self.structure.species: # loop over species
          self.pdos_specie[s] = np.zeros(len(self.energies)) # initialize 
          if self.spin: # spinful calculation
            self.pdos_specie[(s,"up")] = np.zeros(len(self.energies)) 
            self.pdos_specie[(s,"dn")] = np.zeros(len(self.energies)) 
        na = 0 # initialize
        for s in self.structure.species: # loop over species
          for a in self.structure.atomspecies[s]: # loop ove atoms
            na += 1
            self.pdos_specie[s] += self.pdos_atom[a.index] # initialize 
            if self.spin: # spinful calculation
              self.pdos_specie[(s,"up")] += self.pdos_atom[(a.index,"up")] 
              self.pdos_specie[(s,"dn")] += self.pdos_atom[(a.index,"dn")] 
        if na!=len(self.atoms):
            print("Some atom missing in PDOS")
            raise
        # contributions for a specific spin channel
        if self.spin: # spinful calculation
            self.pdos_spin = dict() # dictionary
    def add_total(self,y,spin=None):
        if spin is None: self.dos = y
        elif spin is "up": self.dos_up = y
        elif spin is "dn": self.dos_dn = y
        else: raise

