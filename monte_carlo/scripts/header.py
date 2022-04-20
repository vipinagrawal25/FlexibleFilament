import numpy as np
import matplotlib.pyplot as plt
import sys
import os
Home=os.getenv('Home')
ES=os.getenv('ES')
sys.path.append(ES+'/monte_carlo/tools/')
import functiontools as ft
plt.style.use(Home+'/.matplotlibrc')
plt.rcParams['text.usetex'] = True
import glob
import paramio as pio