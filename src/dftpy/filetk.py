import os
from os import listdir

from os.path import isfile, join


def get_files(mypath,pattern=""):
    """Return files that have a certain string"""
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    out = []
    for f in onlyfiles:
        if pattern in f: out.append(f)
    return out


