# ElasticString
How to run the code:

1> first set an environment variable : ES = <wherver you checked-out the code> 
   This should be done either in .cshrc or .bashrc depending on which shell you are using. 
2> You need to have a file called "host" in the ES directory. This file contains the 
    actual compilation commands that you will use. There are several examples in the
    ES/hosts directory that you can copy or link too. For example, if you are in norlx65
    you can link the file named "hosts/norlx65" to  "host" in the ES directory.
3> Then select a run directory. 
4> At the run directory run:  source $ES/bin/setup.csh
(this is for cshrc, if you use bash you have to write the corresponding bash script)
5> Then copy the two files to your run directory: $ES/Makefile; $ES/input.h 
6> Select the "model" in Makefile
7> Select the input parameters in input.h
8> make 
9> If it compiled then run: ./ode.exe 
