import numpy as np
# import matplotlib.pyplot as plt
import os
os.chdir("../")
import ttisect as T
os.chdir("tests/")
#-------------------------------#
# define triangles
# case -1
# t1 = np.asarray([[0,0,0], [0,1,0], [1,0,0]], dtype=np.float64)
# t2 = np.asarray([[1,0,0], [0,1,0], [1,1,0]], dtype=np.float64)

# case - 2
# t1 = np.asarray([[0,0,0], [0,1,0], [1,0,0]], dtype=np.float64)
# t2 = np.asarray([[0,0,0], [0,1,0], [-1,0,0]], dtype=np.float64)
def read_triangles(file,i):
	t1=np.zeros([3,3])
	t2=np.zeros([3,3])
	t1[0]=file[6*i]
	t1[1]=file[6*i+1]
	t1[2]=file[6*i+2]
	t2[0]=file[6*i+3]
	t2[1]=file[6*i+4]
	t2[2]=file[6*i+5]
	return t1,t2
#
# t2 = t2 + 0.1
# plot the triangles
# def plot_triangles(fig, ax, t1, clr):
#     ax.plot([t1[0][0], t1[1][0]],[t1[0][1], t1[1][1]], 
#             '-', color=clr)
#     ax.plot([t1[0][0], t1[2][0]],[t1[0][1], t1[2][1]], 
#             '-', color=clr)
#     ax.plot([t1[1][0], t1[2][0]],[t1[1][1], t1[2][1]], 
#             '-', color=clr)
#     # ax.plot(t1[2][:-1], t1[0][:-1], '-', color=clr)

# fig, ax = plt.subplots()
# plot_triangles(fig, ax, t1, 'k')
# plot_triangles(fig, ax, t2, 'b')
# plt.show()
lf1=np.loadtxt("rrisect_tri.txt")
size=int(lf1.shape[0]/6)
ncount=0
for i in range(size):
	t1,t2=read_triangles(lf1,i)
	nres=T.tri_tri_isect(t1[0],t1[2],t1[1],t2[1],t2[0],t2[2])
	ncount=ncount+nres
	if nres==1:
		pass
		# print(t1[0,0],t1[0,1],t1[0,2])
		# print(t1[1,0],t1[1,1],t1[1,2])
		# print(t1[2,0],t1[2,1],t1[2,2])
		# #
		# print(t2[0,0],t2[0,1],t2[0,2])
		# print(t2[1,0],t2[1,1],t2[1,2])
		# print(t2[2,0],t2[2,1],t2[2,2])
print(ncount)