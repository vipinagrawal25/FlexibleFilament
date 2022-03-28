import psutil
import sys
sys.path.append("tools/")
import functiontools as ft
import paramio as pio
import numpy as np
from threading import Timer
import time
import os
######################## PARAMETER ###############################
timedelay=10
sigma=0.02
print("# PID = ",os.getpid())
# exit()
######################## SCRIPT ##################################
os.system("make")
# #--------- First stick it ----------------- #
pio.change_param(fname='para_file.in',tip_pos_z=200,is_restart=0,
				afm_sigma=0,afm_epsilon=0,mc_total_iters=10000,mc_dump_iter=10,
				sigma=1.0)
print("# Sticking it")
os.system("mkdir stick")
os.system("./run para_file.in stick > stick/terminal.txt")
ft.wait()
#--------- start simulation for no afm tip ----------------- #
pio.change_param(is_restart=1,sigma=sigma,
				pos_bot_wall=np.loadtxt("stick/mc_log")[-1,-1]-2**(1/6)*sigma)
os.system("mkdir noafm")
os.system("cp stick/restart.h5 noafm/")
print("# No AFM")
os.system("./run para_file.in noafm > noafm/terminal.txt")
ft.wait()
pio.change_param(is_restart=1,afm_sigma=0.02,afm_epsilon=4.0)
#------------------- Now run for all tz --------------------------- #
tz_start=np.loadtxt("noafm/mc_log")[-1,-2]
print("#tz_start\ttz_end")
tz_end=(np.loadtxt("noafm/mc_log")[-1,-1]+tz_start)/2
print(tz_start,"\t",tz_end)
ft.movetip(tz_start,tz_end,restart="noafm")