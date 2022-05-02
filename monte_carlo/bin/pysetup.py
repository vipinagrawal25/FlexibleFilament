import os
import sys
#
ES=sys.path[0]
ES=ES.replace("/bin","/")
#
os.system("mkdir tools")
os.system("mkdir scripts")
os.system("ln -s "+ES+"tools/*.py tools/")
os.system("cp "+ES+"scripts/* scripts/")