compile= nvcc
debug = -dc -g -G 
# options = $(debug)
options = -dc -O3 
link = -l
default: main.o evolve.o model.o cuda.o model.h evolve.h
	$(compile) -o chain.exe evolve.o model.o cuda.o main.o
main.o: main.cu chain.h evolve.h
	$(compile) $(options) main.cu
evolve.o: evolve.cu evolve.h 
	$(compile) $(options) evolve.cu
model.o: model.cu model.h 3vec.h 2Tens.h cuda.h 
	$(compile) $(options) model.cu
cuda.o: cuda.h cuda.cu
	$(compile) $(options) cuda.cu
#------------------------------------------#
tsub: test_sub.o ev.o sub.h   
	$(compile) test_sub.o ev.o 
test_sub.o: test_sub.cu sub.h 
	$(compile) $(options) test_sub.cu 
ev.o: ev.cu sub.h 
	$(compile) $(options) ev.cu 
tall: test_all.o  evolve.o model.o cuda.o
	$(compile) test_all.o evolve.o model.o cuda.o
test_all.o: test_all.cu cuda.h model.h evolve.h 
	$(compile) $(options) test_all.cu 
ttime: test_time.o  evolve.o model.o cuda.o
	$(compile) test_time.o evolve.o model.o cuda.o
test_time.o: test_time.cu  
	$(compile) $(options) test_time.cu 
td: test_diag.o  model.o 
	$(compile) test_diag.o model.o 
test_diag.o: test_diag.cu  
	$(compile) $(options) test_diag.cu 
tp: test_param.cu model.h model.o 
	$(compile)  test_param.cu model.o 
test: test_vec.cu 3vec.h 
	$(compile)  test_vec.cu 
test4.o: test4.cu chain.h
	$(compile) $(options) test4.cu
chain.o: chain.cu chain.h param.h
	$(compile) $(options) chain.cu
tmax: cuda.o test_max.o
	$(compile) test_max.o cuda.o
test_max.o: test_max.cu  
	$(compile) $(options) test_max.cu 
#----------------------------------#
clean:
	rm -f *.o *.exe a.out 
cleann:
	rm -f *.o *.exe a.out core *~ data/*
