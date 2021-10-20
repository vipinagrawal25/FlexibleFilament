import numpy as NP
import os

ShearRate=2;
sigma = [1.5];
sigma = NP.array(sigma);
facAA = [10];
facAA = NP.array(facAA);
dircount = 1;
TMAX = [4];
for ip in range(len(sigma)):
	for jp in range(len(facAA)):
		dirname = "run"+str(dircount)+"/"
		os.system("mkdir " + dirname);
		fparam = open(dirname+"rparam.txt","w+")
		fevolve = open(dirname+"revolve.txt","w+")
		dircount = dircount+1
		fparam.write("%3.2f\t" % sigma[ip])
		fparam.write("%3.2f\t" % facAA[jp])
		fparam.close()
		fevolve.write("%3.2f\t" % TMAX[jp])
		fevolve.write("%d\t" % (TMAX[jp]*50))
		fevolve.close()