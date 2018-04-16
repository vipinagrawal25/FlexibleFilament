mkdir -p src
ln -s $ES/src/Makefile src/
ln -s $ES/src/*.cpp  src/
ln -s $ES/src/*.h src
mkdir -p models
ln -s $ES/src/models/*.cpp models/
ln -s $ES/src/models/*.h models/ 
