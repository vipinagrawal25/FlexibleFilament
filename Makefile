# Makefile for compiling the ode-solver in cpp 
# ---------------------------------------------------------------------
MODEL= estring
#MODEL=particle_in_mag
CUDA=nocuda
target=output
# HOST=norlx65
# ---------------------------------------------------------------------
default:
	(cd src/;ln -sf models/${MODEL}.h model.h;  ln -sf models/${MODEL}.cpp model.cpp; ln -sf ../input.h .; make; mv ode.exe ..)
clean:
	(cd src/; rm -f *.o *.mod *.exe; rm -f model.cpp; rm -f model.h; rm -f input.h;)
multirun:
	(cd src/;ln -sf ../${target}/${MODEL}.h model.h;ln -sf ../${target}/input.h input.h;  ln -sf models/${MODEL}.cpp model.cpp; \
	 make; mv ode.exe ../${target}; cd ..; echo "The executable file has been moved to ${target}.")

