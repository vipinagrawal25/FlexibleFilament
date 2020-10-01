# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
MODEL= estring2dim
#MODEL=particle_in_mag
CUDA=nocuda
folder=output
HOST=norlx65
exec=ode
# ---------------------------------------------------------------------
default:
	(ln -sf hosts/${HOST} host; cd src/;ln -sf models/${MODEL}.h model.h;  ln -sf models/${MODEL}.cpp model.cpp;\
	 ln -sf ../input.h .; make exec=${exec}; mv ${exec}.exe ..)
clean:
	(cd src/; rm -f *.o *.mod *.exe; rm -f model.cpp; rm -f model.h; rm -f input.h;)
multirun:
	(cd src/;ln -sf ../${folder}/${MODEL}.h model.h;ln -sf ../${folder}/input.h input.h;\
	 ln -sf models/${MODEL}.cpp model.cpp; ln -sf ../src/hosts/${HOST} host; make exec=${exec};\
	 mv ${exec}.exe ../${folder}; cd ..; echo "The executable file has been moved to ${folder}.")
dynamics:
	(echo "I will compile myself to compute fixed points or periodic orbits for given bounday condition.\n";\
	 ln -sf hosts/${HOST} host; cd src/;ln -sf models/${MODEL}.h model.h;ln -sf models/${MODEL}.cpp model.cpp;\
	 make map_dyn exec=${exec}; mv ${exec}.exe ..; echo "Compilation completed.")
dynamics_multirun:
	(echo "I will compile myself to compute fixed points or periodic orbits for given bounday condition.\n";\
	 ln -sf hosts/${HOST} host; cd src/;ln -sf ../${folder}/${MODEL}.h model.h;ln -sf models/${MODEL}.cpp\
	 model.cpp; make map_dyn exec=${exec}; mv ${exec}.exe ../${folder}; echo "Compilation completed.")	