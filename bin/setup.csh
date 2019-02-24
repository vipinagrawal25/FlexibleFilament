mkdir -p src
ln -s $ES/src/Makefile src/
ln -s $ES/src/*.cpp  src/
ln -s $ES/src/*.h src
mkdir -p src/models
ln -s $ES/src/models/*.cpp src/models/
ln -s $ES/src/models/*.h src/models/ 
mkdir -p src/modules
ln -s $ES/src/modules/*.h src/modules/ 
cp $ES/Makefile .
cp $ES/input.h .
mkdir -p PythonCodes
ln -s $ES/PythonCodes/*.py PythonCodes/
ln -s $ES/PythonCodes/*.ipynb PythonCodes/
mkdir -p hosts
ln -s $ES/hosts/* hosts/
echo "copied Makefile and input.h.."
echo " Now modify them accroding to your problem."
