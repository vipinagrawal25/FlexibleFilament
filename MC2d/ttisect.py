from ctypes import *
import numpy as np
import os
fpath=os.getcwd()
so_file=fpath+"/include/_libttisect.so"
if os.path.isfile(so_file):
	pass
else:
	os.chdir("include")
	print("Couldn't find .so file for triangle triangle intersection.\nMaking now.")
	os.system("make -s")
_lib=CDLL(so_file)
_lib.tri_tri_intersection_test_3d.argtypes = (POINTER(c_double),POINTER(c_double),POINTER(c_double),
											  POINTER(c_double),POINTER(c_double),POINTER(c_double),
											  c_int,POINTER(c_double),POINTER(c_double))
_lib.tri_tri_intersection_test_3d.restypes = c_int
#
coplanar=0
src=np.zeros(3,dtype=np.double)
tar=np.zeros(3,dtype=np.double)
source=src.ctypes.data_as(POINTER(c_double))
target=tar.ctypes.data_as(POINTER(c_double))
#
def tri_tri_isect(u0,u1,u2,v0,v1,v2):
	vv0=v0.ctypes.data_as(POINTER(c_double))
	vv1=v1.ctypes.data_as(POINTER(c_double))
	vv2=v2.ctypes.data_as(POINTER(c_double))
	uu0=u0.ctypes.data_as(POINTER(c_double))
	uu1=u1.ctypes.data_as(POINTER(c_double))
	uu2=u2.ctypes.data_as(POINTER(c_double))
	return _lib.tri_tri_intersection_test_3d(uu0, uu1, uu2, vv0, vv1, vv2, coplanar,source,target)