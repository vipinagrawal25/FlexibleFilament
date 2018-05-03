# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
compile=g++ -g  -r16
MODEL= harmonic2
#MODEL=particle_in_mag
CUDA=nocuda
# ---------------------------------------------------------------------
default:
	(cd src/;ln -sf models/${MODEL}.h model.h;  ln -sf models/${MODEL}.cpp model.cpp; ln -sf ../input.h .; make; mv ode.exe ..)
clean:
	cd src/; rm -f *.o *.mod *.exe
