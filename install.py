

#import subprocess
#out = subprocess.check_output("which -a python", shell=True)

import os

def get_output(name):
  os.system("rm -f /tmp/command.txt") # remove
  os.system(name+"  > /tmp/command.txt") # run the command
  lines = open("/tmp/command.txt").read() # read the lines
  lines = lines.split("\n") # split the lines
  return lines
  
def get_command(name):
  lines = get_output("which -a "+name) 
  del lines[-1] # remove the last one
  return lines[-1] # return last one

def get_parent(name):
    """Return parent directory"""
    out = get_command(name)
    out = 'echo "$(dirname -- "'+out+'")"'
    return get_output(out)[0]


# now create the symbolic links to Elk and QE
pwd = os.getcwd() # current directory
os.chdir("src/dftpy/codes") # go to folder
os.system("rm -f elk") # remove symbolic link
os.system("ln -s "+get_parent("elk")+"/../ elk")


