compile= nvcc
options = -dc
link = -l
default: odeN.o main.o evolve.o model.o cuda.o
	$(compile) -o odeN.exe odeN.o evolve.o model.o cuda.o main.o
test: cuda.o test.o 
	$(compile)   cuda.o test.o
main.o: main.cu odeN.h evolve.h
	$(compile) $(options) main.cu
test.o: test.cu odeN.h
	$(compile) $(options) test.cu
odeN.o: odeN.cu odeN.h param.h
	$(compile) $(options) odeN.cu
evolve.o: evolve.cu evolve.h 
	$(compile) $(options) evolve.cu
model.o: model.cu model.h 
	$(compile) $(options) model.cu
cuda.o: cuda.h cuda.cu
	$(compile) $(options) cuda.cu
clean:
	rm -f *.o *.exe a.out 
cleann:
	rm -f *.o *.exe a.out core *~
