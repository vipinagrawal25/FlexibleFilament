mkdir -p src
ln -f $ES/src/Makefile src/
ln -f $ES/src/*.cpp  src/
ln -f $ES/src/*.h src
mkdir -p src/models
ln -f $ES/src/models/*.cpp src/models/
ln -f $ES/src/models/*.h src/models/ 
mkdir -p src/modules
ln -f $ES/src/modules/*.h src/modules/ 
cp $ES/Makefile .
ln -f $ES/input.h .
mkdir -p PythonScripts
ln -f $ES/PythonScripts/*.py PythonScripts/
ln -f $ES/PythonScripts/*.ipynb PythonScripts/
echo "copied Makefile and input.h.."
echo " Now modify them accroding to your problem."
