HOST=norlx65_vipin
include hosts/$(HOST)
# CC = h5cc -O3
#
link = -lm -std=c++11
# 
sources = src/forces_lj.c src/forces_surf.c
sources += src/initialise.c  src/visit_io.c src/hdf5_io.c
object =  obj/forces_lj.o obj/initialise.o obj/forces_surf.o
object += obj/visit_io.o obj/hdf5_io.o obj/main.o
# 
includes += include/global.h include/subroutine.h 
# 
run: $(object)
	$(CC) $(object) $(link) -o run
# 
obj/main.o: main.c $(includes)
	$(CC) -Jobj -c main.c -o obj/main.o
# 
object : $(object)
obj/%.o : src/%.c $(includes)
	@mkdir -p $(@D)
	$(info Compiling $<)
	$(CC) -Iobj -c $< -o $@ 
# 
clean:
	@rm $(object) run
	@echo "all object file and executables removed"
# 
distclean:
	@rm -rf data_0* $(object) run
	@echo "all data cleared"