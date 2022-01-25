from ctypes import *
import numpy as np
import os
fpath=os.getcwd()
so_file=fpath+"/tri_tri_intersect.so"
_lib=CDLL(so_file)
_lib.tri_tri_intersection_test_3d.argtypes = (POINTER(c_double),POINTER(c_double),POINTER(c_double),
											  POINTER(c_double),POINTER(c_double),POINTER(c_double),
											 (c_int),POINTER(c_double),POINTER(c_double))
# _lib.tri_tri_intersection_test_3d.restypes = c_int
#
loop_number=10
pair_number=5000
intersection_ratio=0.1
#
def read_triangles(file,i):
	v0=file[6*i]
	v1=file[6*i+1]
	v2=file[6*i+2]
	u0=file[6*i+3]
	u1=file[6*i+4]
	u2=file[6*i+5]
	return v0,v1,v2,u0,u1,u2
#
u0=np.zeros(3)
u1=np.zeros(3)
u2=np.zeros(3)
v0=np.zeros(3)
v1=np.zeros(3)
v2=np.zeros(3)
#Load the triangles
lf1=np.loadtxt("intersected.txt")
lf2=np.loadtxt("separated.txt")
#
nCount=0;
src=np.zeros(3,dtype=np.double)
tar=np.zeros(3,dtype=np.double)
source=src.ctypes.data_as(POINTER(c_double))
target=tar.ctypes.data_as(POINTER(c_double))
coplaner=0
#
nNumber=int(pair_number*intersection_ratio)
for j in range(loop_number):
	for i in range(pair_number):
		if i<nNumber:
			v0,v1,v2,u0,u1,u2=read_triangles(lf1,i)
		else:
			v0,v1,v2,u0,u1,u2=read_triangles(lf2,i)
		vv0=v0.ctypes.data_as(POINTER(c_double))
		vv1=v1.ctypes.data_as(POINTER(c_double))
		vv2=v2.ctypes.data_as(POINTER(c_double))
		uu0=u0.ctypes.data_as(POINTER(c_double))
		uu1=u1.ctypes.data_as(POINTER(c_double))
		uu2=u2.ctypes.data_as(POINTER(c_double))
		nResult=_lib.tri_tri_intersection_test_3d(uu0, uu1, uu2, vv0, vv1, vv2, coplaner,source,target)
		# print(nResult)
		nCount=nCount+nResult
print("nCount=",nCount)