import os
import sys
#
ES=os.getenv('ES')
ES=ES+"/monte_carlo/"
#
if len(sys.argv)==2:
	if sys.argv[1]=='d' or sys.argv[1]=='debug':
		os.system("ln -s" + ES + "main.c .")
		os.system("cp "+ ES +"Makefile .")
		os.system("mkdir hosts")
		os.system("mkdir src")
		os.system("mkdir include")
		os.system("mkdir tools")
		os.system("mkdir input")
		os.system("mkdir run")
		#
		os.system("ln -s "+ES+"src/*.c src/")
		os.system("ln -s "+ES+"include/*.h include/")
		os.system("ln -s "+ES+"tools/*.py tools/")
		os.system("ln -s "+ES+"hosts/* hosts/")
		os.system("cp "+ES+"input/input.h5 input/")
		os.system("cp "+ES+"para_file.in .")
else:
	os.system("cp "+ES+"run .")
	os.system("cp "+ES+"para_file.in .")
	os.system("mkdir input")
	os.system("cp "+ES+"input/*.h5 input/")
	os.system("mkdir run")