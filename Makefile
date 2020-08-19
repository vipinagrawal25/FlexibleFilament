# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
MODEL= estring2dim
#MODEL=particle_in_mag
CUDA=nocuda
target=output
HOST=norlx65
exename=ode
# ---------------------------------------------------------------------
default:
	(ln -sf hosts/${HOST} host; cd src/;ln -sf models/${MODEL}.h model.h;  ln -sf models/${MODEL}.cpp model.cpp; ln -sf ../input.h .;\
	 make exename=${exename}; mv ${exename}.exe ..)
clean:
	(cd src/; rm -f *.o *.mod *.exe; rm -f model.cpp; rm -f model.h; rm -f input.h;)
multirun:
	(cd src/;ln -sf ../${target}/${MODEL}.h model.h;ln -sf ../${target}/input.h input.h;  ln -sf models/${MODEL}.cpp model.cpp;\
	 ln -sf ../src/hosts/${HOST} host; make exename=${exename}; mv ${exename}.exe ../${target}; cd ..; echo "The executable file has been moved to ${target}.")
dynamics:
	(echo "I will compile myself to compute fixed points or periodic orbits for given bounday condition.\n";\
	 ln -sf hosts/${HOST} host; cd src/;ln -sf models/${MODEL}.h model.h;ln -sf models/${MODEL}.cpp model.cpp;\
	 ln -sf ../input.h .; make map_dyn exename=${exename}; mv ${exename}.exe ..; echo "Compilation completed.")