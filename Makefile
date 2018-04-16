# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
compile=g++ -g  -r16
#MODEL=harmonic
MODEL=particle_in_mag
CUDA=nocuda
# ---------------------------------------------------------------------
default:
	(cd src/;ln -sf models/${MODEL}.h model.h;  ln -sf CUDA/${CUDA}.h cuda.h; ln -sf ../input.h .; make; mv ode.exe ..)
clean:
	cd src/; rm -f *.o *.mod *.exe
