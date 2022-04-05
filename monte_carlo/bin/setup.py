import os
import sys
#
ES=sys.path[0]
ES=ES.replace("/bin","/")
# ES=os.getenv('ES')
# ES=ES+"/monte_carlo/"
#
if len(sys.argv)==2:
	if sys.argv[1]=='-d' or sys.argv[1]=='-debug'\
		or sys.argv[1]=='debug' or sys.agrv[1]=='d':
		os.system("ln -s " + ES + "main.c .")
		os.system("ln -s " + ES + "main_mpi.c .")
		os.system("cp "+ ES +"Makefile .")
		os.system("mkdir hosts")
		os.system("mkdir src")
		os.system("mkdir include")
		os.system("mkdir tools")
		os.system("mkdir input")
		os.system("mkdir scripts")
		#
		os.system("ln -s "+ES+"src/*.c src/")
		os.system("ln -s "+ES+"src/*.cpp src/")
		os.system("ln -s "+ES+"src/*.h src/")
		os.system("ln -s "+ES+"include/*.h include/")
		os.system("ln "+ES+"tools/*.py tools/")
		os.system("ln -s "+ES+"hosts/* hosts/")
		os.system("cp "+ES+"input/*.h5 input/")
		os.system("cp "+ES+"para_file.in .")
		os.system("cp "+ES+"scripts/* scripts/")
else:
	curr_dir = os.getcwd()
	os.chdir(ES)
	os.system("make")
	os.chdir(curr_dir)
	os.system("cp "+ES+"run .")
	os.system("cp "+ES+"para_file.in .")
	os.system("mkdir input")
	os.system("cp "+ES+"input/*.h5 input/")
	os.system("cp "+ ES +"Makefile .")
	os.system("mkdir tools")
	os.system("ln "+ES+"tools/*.py tools/")
	os.system("mkdir scripts")
	os.system("cp "+ES+"scripts/* scripts/")