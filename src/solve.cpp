#include <iostream>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "ode.h"
#include "model.h"
#include "cuda.h"
using namespace std;
/**************************/
/* ----------------------------------------*/
int main(){
  double y[ndim];
  double time=0.;
  int k;
  int ldiagnos=0;
  iniconf(y);
  //  setup_cuda();
//----------------------------
  int itn=0;
  while(time <= TMAX){
        printf("%f\t%f\t%f\n",time,y[0],y[1]);
	//    if (ldiagnos==0) {Diagnostics(itn,&y[0],time);};
    for(int ibody=0;ibody<Nensemble;ibody++){
      int irb=pdim*ibody;
      //euler(pdim,&y[irb],time,dt);
      //rnkt2(pdim,&y[irb],time,dt);
      rnkt4(pdim,&y[irb],time,dt);
    }
    time=time+dt;
    itn=itn+1;
    ldiagnos=itn%idiag;
  } 
  cout<<"Done, time="<<time-dt<<"\t TMAX="<<TMAX<<"\n";
//----------------------------
}
/* ----------------------------------------*/

