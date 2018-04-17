# ElasticString
How to run the code:

1> first set an environment variable : ES = <wherver you checked-out the code> 
   This should be done either in .cshrc or .bashrc depending on which shell you are using. 
2> Then select a run directory. 
3> At the run directory run:  source $ES/bin/setup.csh
(this is for cshrc, if you use bash you have to write the corresponding bash script)
4> Then copy the two files to your run directory: $ES/Makefile; $ES/input.h 
5> Select the "model" in Makefile
6> Select the input parameters in input.h
7> make 
8> If it compiled then run: ./ode.exe 
 
