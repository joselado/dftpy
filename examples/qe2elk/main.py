import sys

from dftpy import structure

st = structure.read_qe("qe.in") # read Quantum Espresso structure
structure.write_elk(st) # write Elk structure
