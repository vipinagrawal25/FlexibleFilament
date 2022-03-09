import psutil
import functiontools as f
import paramio as pio
import numpy as np
from threading import Timer
import time
import os
######################## PARAMETER ###############################
tz_start=1.25
tz_end=0.75
npoints=26
timedelay=10
print("# PID = ",os.getpid())
######################## FUNCTIONS ###############################
def isrunning(procname='run'):
	for p in psutil.process_iter(['username','name']):
		if p.info['name']==procname:
			running = 1
		else:
			running = 0
	return running
#---------------------------------------------------------------- #
def wait(procname='run',timedelay=10):
	running=1
	while running==1:
		time.sleep(timedelay)
		running=isrunning()
######################## SCRIPT ##################################
tz_all=np.linspace(tz_start,tz_end,npoints)
# If the code is running, then don't do anything else, change parameters
# and rerun the code.
# os.chdir("../")
# print(os.getcwd())
os.system("make")
#--------- First start simulation for no afm tip ----------------- #
pio.change_param(fname='para_file.in',tip_pos_z=200,is_restart=0,
				afm_sigma=0,afm_epsilon=0)
print("# No AFM")
os.system("mkdir noafm")
os.system("./run para_file.in noafm > noafm/terminal.txt")
wait()
pio.change_param(is_restart=1,afm_sigma=0.17,afm_epsilon=4.0)
#------------------- Now run for all tz --------------------------- #
for tz in tz_all:
	pio.change_param(tip_pos_z=tz)
	print("# tz = ", tz)
	g = float("{0:.3f}".format(tz))
	os.system("mkdir "+str(g))
	os.system("./run para_file.in "+str(g)+" > "+str(g)+"/terminal.txt")
	wait()