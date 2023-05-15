mkdir -p src
ln -s $ES/src/Makefile src/
ln -s $ES/src/*.cpp  src/
ln -s $ES/src/*.h src
mkdir -p src/models
ln -s $ES/src/models/*.h src/models/
ln -s $ES/src/models/*.cpp src/models/
mkdir -p src/modules
ln -s $ES/src/modules/*.h src/modules/
mkdir -p hosts
ln -s $ES/hosts/* hosts/
cp $ES/Makefile .
cp $ES/input.h .
cp $ES/input_mpi.h .

echo "copied Makefile, input.h.. and GPU code"
echo " Now modify them accroding to your problem."