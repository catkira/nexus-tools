#!/usr/bin/env python

import sys
import os.path
import subprocess
import platform

script_path = os.path.dirname(os.path.realpath(__file__))

args = ("python", script_path + "/preprocess.py", "--adapters", "/../data/adapters_best_short.fa", "--flexcat_times", "5", "--flexcat_oh", "5")
args = args + tuple(sys.argv[1:])
popen = subprocess.Popen(args)
popen.wait()
if popen.returncode != 0:
	print "error"
	sys.exit()
