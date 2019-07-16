#compile=g++ -g -Wno-deprecated -std=c++0x
#compile=g++ -pg -O3 -Wno-deprecated -Wno-unused-result -std=c++0x
# CFLAGS=-pg -O3
include ../host
default:
	(${compile} -Wno-unused-result -std=c++0x  -o ode.exe ode.cpp solve.cpp model.cpp)
clean:
	rm -f *.o *.exe 
