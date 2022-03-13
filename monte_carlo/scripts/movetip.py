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
######################## SCRIPT ##################################
os.system("make")
#--------- First start simulation for no afm tip ----------------- #
pio.change_param(fname='para_file.in',tip_pos_z=200,is_restart=0,
				afm_sigma=0,afm_epsilon=0,mc_total_iters=5000,mc_dump_iter=10)
print("# No AFM")
os.system("mkdir noafm")
os.system("./run para_file.in noafm > noafm/terminal.txt")
ft.wait()
pio.change_param(is_restart=1,afm_sigma=0.17,afm_epsilon=4.0)
#------------------- Now run for all tz --------------------------- #
ft.movetip(tz_start,tz_end,npoints)