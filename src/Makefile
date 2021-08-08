#compile=g++ -g -Wno-deprecated -std=c++0x
# compile=g++ -pg -Wno-deprecated -Wno-unused-result -std=c++0x
# CFLAGS=-pg -O3
compile=g++ -I ${ES}/lib -O2 -DNDEBUG -Wno-deprecated -Wno-unused-result -std=c++11
exec=ode
default:
	${compile} -o ${exec}.exe IO.cpp misc.cpp ode.cpp solve.cpp model.cpp
clean:
	rm -f *.o *.exe
map_dyn:
	${compile} -o ${exec}.exe IO.cpp misc.cpp ode.cpp newton-krylov.cpp mapdyn.cpp mapsolve.cpp model.cpp
# 	${compile} -o ${exec}.exe IO.cpp misc.cpp newton-krylov.cpp mapdyn.cpp mapsolve.cpp model.cpp
test1:
	${compile} -o ${exec}.exe misc.cpp newton-krylov.cpp model.cpp test1.cpp mapdyn.cpp
curv2y:
	${compile} -o ${exec}.exe misc.cpp model.cpp test.cpp
tracer:
	${compile} -o ${exec}.exe misc.cpp model.cpp test_tracer.cpp