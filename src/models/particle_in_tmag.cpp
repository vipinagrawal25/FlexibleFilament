#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "ode.h"
#include "3vec.h"
#include "particle_in_mag.h"
using namespace std;
/**************************/
/* ----------------------------------------*/
int main(){
  double y[ndim];
  double time=0.;
  int k;
  int ldiagnos=0;
  iniconf(y);
//----------------------------
  int itn=0;
  while(time <= TMAX){
    //    printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,y[6],y[7],y[8],y[9],y[10],y[11],y[12]);
    if (ldiagnos==0) {Diagnostics(itn,&y[0],time);};
    for(int ibody=0;ibody<Nensemble;ibody++){
      int irb=pdim*ibody;
      rnkt4(pdim,&y[irb],time,dt);
      //euler(pdim,&y[irb],time,dt);
    }
    time=time+dt;
    itn=itn+1;
    ldiagnos=itn%idiag;
  } 
  cout<<"Done, time="<<time-dt<<"\t TMAX="<<TMAX<<"\n";
//----------------------------
}
/* ----------------------------------------*/

