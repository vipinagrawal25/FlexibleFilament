import numpy as np
# def write_param(fname='../para_file.in',N=5120,coef_bending=2.5,coef_stretching=25,
# 				coef_vol_expansion=1e6,radius=1.0,pos_bot_wall=-1.05,sigma=0.17,
# 				epsilon=4,Dfac=32,kBT=1,is_restart=1,mc_total_iters=10000,mc_dump_iter=100,
# 				afm_N=0,extent_l=-0.1,extent_r=0.1,extent_t=-0.1,extent_b=0.1,
# 				tip_radius=0.2,tip_pos_z=1.025,afm_sigma=0.17,afm_epsilon=4):
def write_param(fname="../para_file.in",paramdict=None):
	if paramdict is None:
		paramdict={"N":5120,"coef_bending":2.5,"coef_stretching":25,"coef_vol_expansion":1e6,
				"radius":1.0,"pos_bot_wall":-1.05,"sigma":0.17,"epsilon":4,
				"Dfac":4,"kBT":1,"is_restart":1,"mc_total_iters":10000,"mc_dump_iter":10,
				"afm_N":0,"extent_l":-0.1,"extent_r":0.1,"extent_t":-0.1,"extent_b":0.1,
				"tip_radius":0.2,"tip_pos_z":1.05,"afm_sigma":0.17,"afm_epsilon":4}
	with open(fname, "w") as file:
		file.write("## Membrane parameters\n")
		file.write("N\tcoef_bending\tcoef_stretching\tcoef_vol_expansion\n")
		file.write("%d %3.2f %3.2f %3.2f\n" %(paramdict['N'],paramdict['coef_bending'],
					paramdict['coef_stretching'],paramdict['coef_vol_expansion']))
		file.write("radius\tpos_bot_wall\tsigma\tepsilon\n")
		file.write("%3.2f %3.2f %3.2f %3.2f\n" %(paramdict['radius'],paramdict['pos_bot_wall'],
					paramdict['sigma'],paramdict['epsilon']))
		file.write("## Montecarlo parameters\n")
		file.write("Dfac\tkBT is_restart mc_total_iters\tmc_dump_iter\n")
		file.write("%d %3.2f %d %d %d\n" %(paramdict['Dfac'],paramdict['kBT'],paramdict['is_restart'],
					paramdict['mc_total_iters'],paramdict['mc_dump_iter']))
		file.write("## Afm Tip parameters\n")
		file.write("N\textent_l\textent_r\textent_t\textent_b\n")
		file.write("%d %1.1f %1.1f %1.1f %1.1f\n" %(paramdict['afm_N'],paramdict['extent_l'],
					paramdict['extent_r'],paramdict['extent_t'],paramdict['extent_b']))
		file.write("tip_radius\ttip_pos_z\tafm_sigma\tafm_epsilon\n")
		file.write("%2.2f %2.2f %2.2f %2.2f" %(paramdict['tip_radius'],paramdict['tip_pos_z'],
					paramdict['afm_sigma'],paramdict['afm_epsilon']))
#
def read_param(fname='../para_file.in'):
	## Membrane parameters
	N,coef_bending,coef_stretching,coef_vol_expansion = np.loadtxt(fname,skiprows=2,max_rows=1)
	N=int(N)
	radius,pos_bot_wall,sigma,epsilon = np.loadtxt(fname,skiprows=4,max_rows=1)
	## Montecarlo parameters
	Dfac,kBT,is_restart,mc_total_iters,mc_dump_iter = np.loadtxt(fname,skiprows=7,max_rows=1)
	Dfac=int(Dfac)
	is_restart=int(is_restart)
	mc_total_iters=int(mc_total_iters)
	mc_dump_iter=int(mc_dump_iter)
	## Afm Tip parameters
	afm_N,extent_l,extent_r,extent_t,extent_b = np.loadtxt(fname,skiprows=10,max_rows=1)
	afm_N=int(afm_N)
	tip_radius,tip_pos_z,afm_sigma,afm_epsilon = np.loadtxt(fname,skiprows=12,max_rows=1)
	paramdict={"N":N,"coef_bending":coef_bending,"coef_stretching":coef_stretching,
				"coef_vol_expansion":coef_vol_expansion,
				"radius":radius,"pos_bot_wall":pos_bot_wall,"sigma":sigma,"epsilon":epsilon,
				"Dfac":Dfac,"kBT":kBT,"is_restart":is_restart,"mc_total_iters":mc_total_iters,
				"mc_dump_iter":mc_dump_iter,
				"afm_N":afm_N,"extent_l":extent_l,"extent_r":extent_r,"extent_t":extent_t,
				"extent_b":extent_b,
				"tip_radius":tip_radius,"tip_pos_z":tip_pos_z,"afm_sigma":afm_sigma,
				"afm_epsilon":afm_epsilon}
	return paramdict
#
def change_param(fname="../para_file.in",**kwargs):
	''' The function change one parameter from the input file and overwrites the new parameter file.'''
	paramdict=read_param(fname=fname)
	for key,value in kwargs.items():
		paramdict[key]=value
	write_param(fname=fname,paramdict=paramdict)
#