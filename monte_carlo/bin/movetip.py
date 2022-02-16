import psutil
import funcbin as f
import numpy as np
from threading import Timer
import time
import os
########################PARAMETER###############################
tz_start=1.05
tz_end=0.75
npoints=13
timedelay=1
########################FUNCTIONS###############################
def isrunning(pname='run'):
	for p in psutil.process_iter(['username','name']):
		if p.info['name']==pname:
			running = 1
		else:
			running = 0
	return running
########################SCRIPT###############################
tz_all=np.linspace(tz_start,tz_end,npoints)
# If the code is running, then don't do anything else, change parameters
# and rerun the code.
os.chdir("../")
# print(os.getcwd())
for tz in tz_all:
	running=1
	while running==1:
		time.sleep(timedelay)
		running=isrunning()
	f.change_param(fname='para_file.in',tip_pos_z=tz,is_restart=1)
	g = float("{0:.3f}".format(tz))
	os.system("./run para_file.in "+str(g))