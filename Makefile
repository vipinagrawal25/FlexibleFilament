# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
MODEL=shell
#MODEL=particle_in_mag
CUDA=nocuda
folder=output
HOST=norlx65
exec=ode
# ---------------------------------------------------------------------
default:
	(cd src/; ln -sf ../hosts/${HOST} host; ln -sf models/${MODEL}.h model.h; ln -sf models/${MODEL}.cpp model.cpp;\
	 ln -sf ../input.h .; make exec=${exec}; mv ${exec}.exe ..)
clean:
	(cd src/; rm -f *.o *.mod *.exe; rm -f model.cpp; rm -f model.h; rm -f input.h;)
multirun:
	(cd src/;ln -sf ../hosts/${HOST} host; ln -sf ../${folder}/${MODEL}.h model.h;ln -sf ../${folder}/input.h input.h;\
	 ln -sf models/${MODEL}.cpp model.cpp; make exec=${exec};\
	 mv ${exec}.exe ../${folder}; cd ..; echo "The executable file has been moved to ${folder}.")
dynamics:
	(echo "I will compile myself to compute fixed points or periodic orbits for given bounday condition.\n";\
	 cd src/;ln -sf ../hosts/${HOST} host; ln -sf models/${MODEL}.h model.h;ln -sf models/${MODEL}.cpp model.cpp;\
	 make map_dyn exec=${exec}; mv ${exec}.exe ..; echo "Compilation completed.")
dynamics_multirun:
	(echo "I will compile myself to compute fixed points or periodic orbits for given bounday condition.\n";\
	 cd src/; ln -sf ../hosts/${HOST} host; ln -sf ../${folder}/${MODEL}.h model.h;ln -sf models/${MODEL}.cpp\
	 model.cpp; make map_dyn exec=${exec}; mv ${exec}.exe ../${folder}; echo "Compilation completed.")
mpi:
	(cd src/; ln -sf ../hosts/${HOST} host;ln -sf models/${MODEL}_mpi.h model.h;ln -sf models/${MODEL}.cpp model.cpp;\
	 ln -sf ../input_mpi.h input.h; make MPI exec=${exec}; mv ${exec}.exe ..)