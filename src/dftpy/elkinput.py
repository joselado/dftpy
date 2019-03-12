from . import units

def extract_attributes(cal):
    """
    Return a string with attributes that need to be written to elk.in
    """
    out = ""
    out += "epsengy\n" + str(cal.energyerror/units.h2ev)+"\n\n"
    return out


def generate_input(d):
    """
    Add additional input variables to the Elk file based on a dictionary
    """
    out = "" # output
    for di in d: # loop over dictionary
        o = python2elk(d[di])
        out += di+"\n" # add string
        out += o+"\n\n\n" # add string
    return out # return output



def python2elk(a):
    """
    Transform a Python variable in an Elk variable
    """
    if type(a) is type(True): # bool variable
        if a: return ".true."
        else: return ".false."
    else: return str(a)

