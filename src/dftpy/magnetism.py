import numpy as np

datoms = ["Ni","Cu","Fe"]

def lowest_energy_state(c0,atoms=datoms):
    """Return the configuration with lowest energy"""
    cs = [] # empty list
    c = c0.copy() # copy calculation
    c.spin = True # set spinful
    c.struct.set_magnetism(m) # set the magnetic configuration
    e = c.gs_energy() # get the GS energy
    return cmin # return


