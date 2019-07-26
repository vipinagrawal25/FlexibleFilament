compile= nvcc
options = -dc 
link = -l
default: chain.o main.o evolve.o model.o cuda.o
	$(compile) -o chain.exe chain.o evolve.o model.o cuda.o main.o
test: model.o test4.o 
	$(compile)   model.o test4.o
main.o: main.cu chain.h evolve.h
	$(compile) $(options) main.cu
test4.o: test4.cu chain.h
	$(compile) $(options) test4.cu
chain.o: chain.cu chain.h param.h
	$(compile) $(options) chain.cu
evolve.o: evolve.cu evolve.h 
	$(compile) $(options) evolve.cu
model.o: model.cu model.h 3vec.h 2Tens.h 
	$(compile) $(options) model.cu
cuda.o: cuda.h cuda.cu
	$(compile) $(options) cuda.cu
clean:
	rm -f *.o *.exe a.out 
cleann:
	rm -f *.o *.exe a.out core *~
