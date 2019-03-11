#!/usr/bin/python

import os

os.system("python clean.py")
os.system("git add .")
os.system("git commit -m 'Update'")
os.system("git push")

