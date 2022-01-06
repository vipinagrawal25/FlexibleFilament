compile= g++ -std=c++11
debug = -dc -g -G
# options = $(debug)
options = -dc -O3 -w
link = -l
default: 
	$(compile) model.cpp solve.cpp
clean:
	rm -f *.o *.exe a.out 
#----------------------------------#