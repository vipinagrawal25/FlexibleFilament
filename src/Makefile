#compile=g++ -g -Wno-deprecated -std=c++0x
# compile=g++ -pg -Wno-deprecated -Wno-unused-result -std=c++0x
# CFLAGS=-pg -O3
include ../host
exec=ode
default:
	${compile} -std=c++11 -o ${exec}.exe ode.cpp solve.cpp model.cpp
clean:
	rm -f *.o *.exe
map_dyn:
	${compile} -std=c++11 -o ${exec}.exe ode.cpp map_dyn.cpp model.cpp