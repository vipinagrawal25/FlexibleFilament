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

### Configuration Number
This section is about defining conf_number in src.h file

1> Conf_Number == 0 :- A rod which is at equilibrium in vertical direction but have a small deformation in starting in xy plane. This configuration is fixed at bottom.
2> Conf_Number == 1 :- This configuration is used for demonstration of GI Taylor's experiment. In this case, the rod is free and nowhere clamped. Equation of the rod is X = 0. In this case the shear flow is in perpendicular direction to the rod .
3> Conf_Number == 2 :- In this case, a free rod is kept in the direction of the flow (at the origin so no force) but the rod is given a small deviation in starting we want to study the dynamics of that.

> Tip: Define lastfile variable to a non-zero value to start the code from middle.