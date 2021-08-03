#include <iostream>
#include <fstream>
#include "modules/2vec.h"
#include "modules/3vec.h"
#include "modules/2b2Tens.h"
#include "modules/2Tens.h"
#include "model.h"
#include <string>
#include <math.h>
#include "constant.h"
#include "IO.h"
#include <memory.h>
#include "misc.h"
#include "MapDyn.h"
/**************************/
using namespace std;
/**************************/
/* 1) We consider simulation of elastic filament in 2 dimension. Filament is in XY dimension.
      X direction is horizontal. */
/**************************/
void dHdR(int kp, vec2 X[], vec2* add_FF, double* add_kappasqr, bool flag_kappa,double time);
void getub(double *bk, vec2 *uk, int kp, vec2 X[]);
int MatrixtoVector(int i, int j, int N);
void GetRij(vec2 X[], int i, int j, double *Distance, vec2 *rij);
void drag(vec2 X[], vec2 dX[], vec2 EForce[]);
vec2 y2vec2(double *y,int ip);
vec3 y2vec3(double *y,int ip);
void ext_force(vec2* EForce,double *y,double time);
void vec2y(double *y, vec2 XX, int ip);
void vec2y(double *y, vec3 XX, int ip);
vec3 ext_flow(vec3 R3d, double time);
vec2 ext_flow(vec2 R, double time);
void calc_Xconstrain(vec2* Xcons, double time);
/**************************/
// Global variables to the file.
/* save the coordinates of first two points to go back and forth from kappa to y.
   first 4 points for points before map iteration and last 4 points for after map iteration. */
double y_start[2*pp] = {0.,0.,0.,aa};
string cc;
double Error=0;
/**************************/
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
  double EForceArr[ndim];
  eval_rhs(time,y,rhs,flag_kappa,CurvSqr,SS,EForceArr,0);
}
/**************************/
void eval_rhs(double time, double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[], 
              double EForceArr[], bool iEforceArr){
  vec2 R[Np],dR[Np],EForce_ip,FF0,EForce[Np];
  // R is the position of the beads.
  // double CurvSqr[Np];
  double kappasqr;
  double onebythree = 1./3.;
  SS[0] = 0;                                      // Initializing the material co-ordinate
  for (int ip=0;ip<Np;ip++){
    R[ip].x=y[2*ip];
    R[ip].y=y[2*ip+1];
  }
  for (int ip=0;ip<Np;ip++){
    kappasqr=CurvSqr[ip];
    if (flag_kappa){
      if (ip<Np-1){
        SS[ip+1] = SS[ip] + norm(R[ip+1]-R[ip]);
      }
    }
    dHdR(ip, R, &EForce_ip, &kappasqr, flag_kappa, time);
    EForce[ip] = EForce_ip;
    CurvSqr[ip]=kappasqr;
  }
  // Is the string forced at some points? Implement it here.
  
  ext_force(EForce,y,time); 
  drag(R, dR, EForce);
  for (int ip = 0; ip < Np; ++ip){
    dR[ip]=dR[ip]+ext_flow(R[ip],time);
  }
  for (int ip=0;ip<Np;ip++){
    if (iEforceArr){
      EForceArr[2*ip] = EForce[ip].x;
      EForceArr[2*ip+1] = EForce[ip].y;
    }
    rhs[2*ip]=dR[ip].x;
    rhs[2*ip+1]=dR[ip].y;
  }
}
/**************************/
void eval_rhs_tr(double time, double EForceArr[], double y[], double y_tr[], double rhs[]){
  int ndim_tr=np_tracer*pp_tracer;
  vec3 dX, RR, EForce[Np], Rtracer;
  vec2 FF(0.,0.), R[Np];
  for (int ip = 0; ip < Np; ++ip){
    EForce[ip].x = 0;
    EForce[ip].y = EForceArr[2*ip]-FF.x;
    EForce[ip].z = EForceArr[2*ip+1]-FF.y;
  }
  double d_RR,c1,dsqr1;
  Tens2 mu;
  double onebythree = 1./3.;
  for (int ip_tracer = 0; ip_tracer < np_tracer; ++ip_tracer){
    dX.x = 0;
    dX.y = 0;
    dX.z = 0;
    Rtracer = y2vec3(y_tr, ip_tracer);
    for (int ip = 0; ip < Np; ++ip){
      RR.x = Rtracer.x;
      RR.y = Rtracer.y - y[2*ip];
      RR.z = Rtracer.z - y[2*ip+1];
      d_RR = norm(RR);
      c1 = 1/(8*M_PI*viscosity*d_RR);
      dsqr1 = 1./(d_RR*d_RR);
      mu = c1*(dab + (RR*RR)*dsqr1 + 1./4.*dd*dd*dsqr1*(dab*onebythree - (RR*RR)*dsqr1));
      dX = dX + dot(mu,EForce[ip]);
      // PVec3(dX);
    }
    // if we want to save the velocities and positions without external shear flow.
    // vec2y(rhs, dX, ip_tracer+np_tracer);
    dX = dX + ext_flow(Rtracer, time);
    vec2y(rhs, dX, ip_tracer);
  }
}
/**************************/
void ext_force(vec2* EForce, double* y, double time){
  vec2 Xghost[np_cons],XX;
  double GG = 10000;
  if (icons){
    calc_Xconstrain(Xghost,time);
    // PVec2(Xghost[0]);
    for (int i = 0; i < np_cons; ++i){
      XX = y2vec2(y,loc_con[i]);
      EForce[loc_con[i]] = EForce[loc_con[i]] - (XX-Xghost[i])*GG;
    }
  }
}
/**************************/
vec2 ext_flow(vec2 R_ip, double time){
  double omega = 2*M_PI/period;
  vec2 dR;
  double Ts=40;
  switch(iext_flow){
    case 1:
    if (sin(omega*time)>=0){
        dR.x = dR.x + ShearRate*(height-R_ip.y)*ceil(sin(omega*time));
    }
    else{
        dR.x = dR.x + ShearRate*(height-R_ip.y)*floor(sin(omega*time));
    }
    break;
    case 2:
      // NOT SO GOOD WAY TO WRITE LIKE THIS
      if (time<Ts*(1/ShearRate)){
        dR.x = dR.x + ShearRate*(R_ip.y);
      }else{
        dR.x = dR.x + ShearRate*(R_ip.y)*sin(omega*(time-Ts*(1/ShearRate)));
      }
      break; 
    case 3:
      dR.x = dR.x + ShearRate*(height-R_ip.y)*sin(omega*time);
      break;
  }
  return dR;
}
/**************************/
vec3 ext_flow(vec3 R3d, double time){
  // This is a wrapper function for vec3 inputs.
  vec2 R2d;
  vec3 dR;
  R2d.x = R3d.y;
  R2d.y = R3d.z;
  //
  R2d = ext_flow(R2d,time);
  dR.x = 0;
  dR.y = R2d.x;
  dR.z = R2d.y;
  return dR; 
}
/**************************/
void drag(vec2 X[], vec2 dX[], vec2 EForce[]){
  double onebythree = 1./3.;
  double mu0 = onebythree/(M_PI*viscosity*dd);
  if (UseRP == 'Y'){
    // mu_ij represents the one element of mobility matrix (Size: NXN).
    // Every element of the matrix itself is a 2nd rank tensor and the dimension of that should 3x3.
    Tens2b2 mu_ij, mu_ii;
    double d_rij;
    vec2 rij;
    mu_ii = dab2b2*mu0;       // dab``2b2 is the unit 2 dimensional tensor. It is defined in module/2b2Tens file.
    for (int ip = 0; ip < Np; ++ip){
      // The mu_ij in the next line represents the mobility tensor when j is equal to i and in 
      // response to that the "for loop" is started from ip+1.
      dX[ip] = dX[ip] + dot(mu_ii, EForce[ip]);
      for (int jp = ip+1; jp < Np; ++jp){
        GetRij(X, ip, jp, &d_rij, &rij);
        double c1 = 1/(8*M_PI*viscosity*d_rij);
        double dsqr1 = 1./(d_rij*d_rij);
        mu_ij = c1*(dab2b2 + (rij*rij)*dsqr1 + 1./2.*dd*dd*dsqr1*(dab2b2*onebythree - (rij*rij)*dsqr1));
        dX[ip] = dX[ip] + dot(mu_ij, EForce[jp]);
        dX[jp] = dX[jp] + dot(mu_ij, EForce[ip]);
      }
    }
  }
  else if(UseRP=='N'){
    for (int ip = 0; ip < Np; ++ip){
        dX[ip] = EForce[ip]*mu0;
    }
  }
  else{
    cout << "Sorry, I do not understand your choice of the paramenter UseRP." << endl;
    exit(1);
  }
}
/**************************/
void getub(double *bk, vec2 *uk, int kp, vec2 X[]){
  vec2 dX = X[kp+1]-X[kp];
  double bb = norm(dX);
  *bk = bb;
  *uk =dX/bb;
}
/**************************/
vec2 extension_force(int kp, vec2 X[], double time){
  vec2 FF;
  vec2 ukm1(0.,0.), uk(0.,0.), Xzero(0.,0.), dX(0.,0.), Xone(0.,0.);
  double yzero[2] = {0.,0.};
  double yone[2] = {0.,0.};
  double bkm1,bk;
  if (bcb==0){
    // If the bottom point is fixed, then it can not move.
    Xzero.x=0.; Xzero.y=0.;
  }else if(bcb==2){
    // We define the coordinates for clamped boundary condition also.
    Xzero.x=0; Xzero.y=0;
    calc_yone(yone,time);
    Xone.x = yone[0]; Xone.y = yone[1];
  }
  switch(kp){
    case 0:
      getub(&bk, &uk, kp, X);
      if (bcb==0){
          dX = X[kp-1+1]-Xzero;
          bkm1 = norm(dX);
          ukm1=dX/bkm1;
      }
      else if(bcb==1){}
      else if(bcb==2){
        dX = X[0]-Xone;
        bkm1 = norm(dX);
        ukm1=dX/bkm1;
      }
      else{
        cout << "Boundary condition at bottom not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      break;

    case Np-1:
      if(bct==1){
        getub(&bkm1, &ukm1, kp-1, X);
      }
      else{
        cout << "Boundary condition at Np==NN not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      break;

    default:
      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);      
      break;
  }
  // FF = extension_force(ukm1,uk,bkm1,bk);
  FF = (uk*(bk-aa) - ukm1*(bkm1-aa))*HH/aa;
  return FF; 
}
/**************************/
// vec2 extension_force(vec2 ukm1, vec2 uk, double bkm1, double bk){
//   vec2 FF;
//   FF = (uk*(bk-aa) - ukm1*(bkm1-aa))*HH/aa;
//   return FF;
// }
/**************************/
// vec2 bending_force(vec2 ukm2, vec2 ukm1, vec2 uk, vec2 ukp1, double bkm2, double bkm1, double bk, double bkp1){
//   vec2 FF;
//   if(bk<tiny){
//     FF =  ukm2/bkm1- (ukm1/bkm1)*( dot(ukm1,uk) ) ;
//   }
//   else if(bkm1<tiny){
//     FF =  (uk/bk)* dot(uk,ukp1)  - (ukm1+ukp1)/bk ;
//   }
//   else{
//     FF = ( (uk+ukm2)/(bkm1+tiny) - (ukm1+ukp1)/(bk+tiny)
//        + (uk/(bk+tiny))*( dot(uk,ukm1) + dot(uk,ukp1) )
//        - (ukm1/(bkm1+tiny))*( dot(ukm1,uk) + dot(ukm1,ukm2) )
//       );
//   }
//   FF = FF*AA/aa;
//   return FF;
// }
/***********************/
double square_wave(double time,double period){
  int s;
  double sw;
  s = floor(2*time/period);
  sw = (double)-2*(s%2)+1;
  return sw;
}
/**************************/
void calc_Xconstrain(vec2* Xcons, double time){
  double inidis = 0;
  // double period = 2*M_PI/omega;
  switch(icons){
    case 1:
      // points on circle
      for(int i = 0; i < np_cons; ++i){
        Xcons[i].x = (inidis +  aa/height*(loc_con[i]-loc_con[0]))*cos(angularVel*time);
        Xcons[i].y = (inidis +  aa/height*(loc_con[i]-loc_con[0]))*sin(angularVel*time);
      }
      break;
    case 2:
      // points on figure '8'
      for (int i = 0; i < np_cons; ++i){
        Xcons[i].x = (height +  aa*(loc_con[i]-loc_con[0]))*cos(angularVel*time);
        Xcons[i].y = (height +  aa*(loc_con[i]-loc_con[0]))*cos(angularVel*time)*sin(angularVel*time);
      }
    case 3:
      // angularVel*time should change its direction.
      // double theta = 4*M_PI*abs(time/period - floor(time/period+0.5));
      double theta = angularVel*time*square_wave(time,period);
      for (int i = 0; i < np_cons; ++i){
        Xcons[i].x = (inidis +  aa/height*(loc_con[i]-loc_con[0]))*cos(theta);
        Xcons[i].y = (inidis +  aa/height*(loc_con[i]-loc_con[0]))*sin(theta);
      }
      break;
    case 4:
      for (int i = 0; i < np_cons; ++i){
        Xcons[i].x = (inidis +  aa/height*(loc_con[i]-loc_con[0]))*cos(theta);
        Xcons[i].y = (inidis +  aa/height*(loc_con[i]-loc_con[0]))*sin(theta);
      }
    }
}
/**************************/
void calc_yone(double *yone, double time){
  double angularVel = 2*M_PI/period;
  double yzero[2];
  calc_yzero(yzero,time);
  yone[0] = yzero[0] + aa*cos(angularVel*time);
  yone[1] = yzero[1] + aa*sin(angularVel*time);
  // print(yone,2);
}
/**************************/
void calc_yzero(double *yzero, double time){
  double angularVel = 2*M_PI/period;
  yzero[0] = height/2*cos(angularVel*time);
  yzero[1] = height/2*sin(angularVel*time);
}
/**************************/
void dHdR(int kp, vec2 X[], vec2* add_FF, double* add_kappasqr, bool flag_kappa, double time){
  // This function calculates the force at every node which is a function of X, time.
  vec2 ukm2(0.,0.), ukm1(0.,0.), uk(0.,0.), ukp1(0.,0.), Xzero(0.,0.), dX(0.,0.), Xone(0.,0.);
  double bkm2 = 0.;
  double bkm1 = 0.;
  double bk = 0.;
  double bkp1 = 0.;
  vec2 FF = *add_FF;
  double yzero[2] = {0.,0.};
  double yone[2] = {0.,0.};
  // Since I am passing the address of force in add_FF and the same goes for Kapppsqr
  //vec2 FF;
  /* Here the problem is that Xzero has been taken as the first point of the rod and which is claimed to be 
  fixed in general. But for some cases like the implementation of the taylor experiment, 
  we want this to be free. For this I am implementing Xzero based on the configuration. 
  Since finally with this function we just want to calculate the force on particular node. 
  So we can just change the way we calculate force on 1st and 2nd node for different configuration.*/
  /* One thing should be noted that now the expression inside case 0 and case 1 would change depending upon 
  the configuration. Suppose if XZero is taken as the fixed point, then we dont need to calculate F^{0} 
  but only F^{1} which would be implemented in case 0 because Xzero is the bottom most point for 
  which we dont care to calculate the force and X[0] is the first point being implemented in case 0. 
  Though for a different configuration these things should be changed. */
  switch(kp){
    case 0:
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      if (bcb==0){
        Xzero.x=0.; Xzero.y=0.;
        dX = X[kp-1+1]-Xzero;
        bkm1 = norm(dX);
        ukm1=dX/bkm1;
        //
        FF = (  (uk)/bkm1 - (ukm1+ukp1)/bk
             + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
             - (ukm1/bkm1)*( dot(ukm1,uk) )
             );
        FF = FF*AA/aa;
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;
      }
      else if(bcb==1){
        FF =  (uk/bk)*dot(uk,ukp1)  - (ukp1)/bk ;
        FF = FF*AA/aa;
        // Add an extra term for inextensibility constraint
        FF = FF + (uk*(bk-aa))*HH/aa;
      }
      else if (bcb==2){
        calc_yzero(yzero,time);
        calc_yone(yone,time);
        Xone.x = yone[0]; Xone.y = yone[1];
        Xzero.x = yzero[0]; Xzero.y = yzero[1];
        //
        dX = X[0]-Xone;
        bkm1 = norm(dX);
        ukm1=dX/bkm1;
        //
        dX = Xone-Xzero;
        bkm2 = norm(dX);
        ukm2=dX/bkm2;
        //
        FF = ( (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
        FF = FF*(AA/aa);
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;            // Inextensibility constraint
      }
      else{
        cout << "Boundary condition at bottom not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      *add_kappasqr=0.;
      break;
    case 1:
      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      if (bcb==0){
        Xzero.x=0.; Xzero.y=0.;
        dX = X[kp-2+1]-Xzero;
        bkm2 = norm(dX);
        ukm2 = dX/bkm2;
        //
        FF = ( (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
        FF = FF*(AA/aa);
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;        // Inextensibility constraint
      }else if(bcb==1){
        FF = ( (uk)/bkm1 - (ukm1+ukp1)/bk
              + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
              - (ukm1/bkm1)*( dot(ukm1,uk) )
              );
        FF = FF*(AA/aa);
        // cout << FF.z << endl;
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;        // Inextensibility constraint
      }
      else if(bcb==2){
        calc_yone(yone,time);
        Xone.x = yone[0]; Xone.y = yone[1];
        //
        dX = X[0]-Xone;
        bkm2 = norm(dX);
        ukm2=dX/bkm2;
        //
        FF = ( (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
        FF = FF*(AA/aa);
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;          // Inextensibility constraint
      }
      else{
        cout << "Boundary condition at Np==1 not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      *add_kappasqr=0.;
      break;
    case Np-2:
      if(bct==1){
        getub(&bkm2, &ukm2, kp-2, X);
        getub(&bkm1, &ukm1, kp-1, X);
        getub(&bk, &uk, kp, X);
        *add_kappasqr=0.;
        FF = (     (uk+ukm2)/bkm1 - (ukm1)/bk
            + (uk/bk)*( dot(uk,ukm1))
            - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
            );    
        FF = FF*(AA/aa);
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
      }
      else{
        cout << "Boundary condition at Np==NN-1 not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      break;
    case Np-1:
      if(bct==1){
        getub(&bkm2, &ukm2, kp-2, X);
        getub(&bkm1, &ukm1, kp-1, X);
        *add_kappasqr=0.;
        //
        FF = (ukm2)/bkm1 - (ukm1/bkm1)*dot(ukm1,ukm2) ;
        FF = FF*(AA/aa);
        FF = FF - (ukm1*(bkm1-aa))*HH/aa;
      }
      else{
        cout << "Boundary condition at Np==NN not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      break;
    default:
      getub(&bkm2, &ukm2, kp-2, X);
      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      if (flag_kappa){
        *add_kappasqr=2.*(1.- dot(uk,ukm1))/(aa*aa);
      }
      FF = ( (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
      FF = FF*(AA/aa);
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;          // Inextensibility constraint
      break;
  }
  *add_FF = FF;
}
/****************************************************/
void iniconf(double *y, double *aTime, double tdiag){
  int itn = (int) ((*aTime)/tdiag);
  cout << itn << endl;
  string fname;
  if (FileExists("output/var"+ to_string(itn) + ".txt")){
    fname = "output/var" + to_string(itn) + ".txt";
  }
  else if(FileExists("var" + to_string(itn) + ".txt")){
    fname = "var"+ to_string(itn) + ".txt";
  }
  else{
    cout << "# Input file does not exist for time = " << *aTime << endl;
    cout << "# I shall start the simulation from time = 0 " << endl;
    *aTime = 0;
    exit(1);
  }
  rData(y,fname);
}
/****************************************************/
void iniconf(double *y){
  // Merge all three cases into 1;
  vec2 R[Np];              // R is the position of the beads.
  double k = 2;            // determines the frequency for initial configuration
  double CurvLength = 0;   // determines the total length of the curve
  string l;
  // double theta= (double)M_PI/2;
  double theta = 0;
  // double vdis=2*aa;      // Define a parameter called middle point in model.h
  double inidis=0;          // that will take care of everything
  //
  vec2 Xcons[np_cons],XX;
  ifstream myfile;
  double ch;
  int cnt=0;
  switch(niniconf){
    case -2:
      // To give perturbation to the filament.
      rData(y,datafile);
      print(y,Np,pp);
      cout << endl;
      // add perturbation to the file.
      for (int ip = 0; ip < Np; ++ip){
        XX = y2vec2(y,ip);
        XX.x = XX.x + 5*aa*sin(M_PI*k*aa*double(ip)/height);
        vec2y(y, XX, ip);
      }
      print(y,Np,pp);
      break;
    case -1:
      // myfile.open(datafile);
      // // myfile >> ch;
      // while(myfile >> ch){cnt++;}
      // myfile.close();
      // cout << "# Number of elements read: " << cnt << endl;
      // if(cnt==Np){}
      // else if(cnt==ndim){
      //   cnt=0;
      //   myfile.open(datafile);
      //   while(myfile >> ch){y[cnt]=ch;cnt++;}
      //   memcpy(y_start,y,2*pp*sizeof(double));
      //   myfile.close();
      // }
      // else if(cnt==pp*ndim){
      //   cout << datafile[0] << endl;
      //   myfile.open(datafile);
      //   for(int idim = 0; idim < Np; ++idim){
      //     for (int ip = 0; ip < pp; ++ip){
      //       myfile >> y[idim*pp+ip];
      //     }
      //     for (int ip = 0; ip < pp; ++ip){
      //       myfile >> ch;
      //     }
      //   }
      //   memcpy(y_start,y,2*pp*sizeof(double));
      //   myfile.close();
      // }
      // else{
      //   cout << "# Initial point is not right." << endl;
      //   exit(1);
      // }
      // myfile.close();
      rData(y,datafile);
      break;
    case 0:
      y[0]=0;
      y[1]=0;
      for (int ip=1;ip<Np;ip++){
        R[ip].x=aa*sin(M_PI*k*aa*double(ip)/height);
        R[ip].y=aa*double(ip);
        CurvLength += norm(R[ip]-R[ip-1]);
        y[2*ip] = R[ip].x;
        y[2*ip+1] = R[ip].y;
      }
      for (int idim = 0; idim < ndim; ++idim){
        y[idim] = y[idim]*height/CurvLength;
      }
      break;
    case 1:
      // In this case we implement the initial configuration for GI Taylor experiment. 
      // i.e. a straight rod which has length equal to the height of the box and free to move from bottom.
      for (int ip = 0; ip < Np; ++ip){
        R[ip].x = (aa*double(ip)+inidis)*sin(theta);
        R[ip].y = (aa*double(ip)+inidis)*cos(theta);
        //
        y[2*ip] = R[ip].x;
        y[2*ip+1] = R[ip].y;
      }
      break;
    case 2:
      // Santillian's experiment
      // In this case, we want to study the dynamics of a rod which is kept in the direction of the flow at 
      // origin.  The rod should be deviated a little bit from origin in starting.
      for (int ip = 0; ip < Np; ++ip){
        R[ip].x = aa*(double(ip+1))-height*0.5;
        R[ip].y = aa*sin(M_PI*k*aa*double(ip+1)/height);
        if (ip>0){
          CurvLength = CurvLength + norm(R[ip] - R[ip-1]);
            // cout << CurvLength << endl;
        }
        else{
          CurvLength = CurvLength + sqrt((R[ip].x)*(R[ip].x)+(R[ip].y)*(R[ip].y));
        }
        y[2*ip] = R[ip].x;
        y[2*ip+1] = R[ip].y;
      }
      for (int ip = 0; ip < Np; ++ip){
        y[2*ip+1] = R[ip].y/CurvLength;
      }
      break;
    case 3:
      calc_Xconstrain(Xcons,0);
      for (int ip_cons = 0; ip_cons < np_cons; ++ip_cons){
        R[loc_con[ip_cons]] = Xcons[ip_cons];
      }
      if(norm(R[loc_con[np_cons-1]] - R[loc_con[0]]) < aa*loc_con[np_cons-1] - aa*loc_con[0]){
         k = height/(aa*loc_con[np_cons-1] - aa*loc_con[0]);
      }
      //
      for (int ip = 0; ip < Np; ++ip){
        R[ip].x = double(ip)/(loc_con[np_cons-1]-loc_con[0]);
        R[ip].y = sin(M_PI*k*R[ip].x);
        if (ip>0){
          CurvLength = CurvLength + norm(R[ip] - R[ip-1]);
        }
      }
      //
      while(abs(CurvLength-height) > 0.01){
        for (int ip = 0; ip < Np; ++ip){
          R[ip].y = R[ip].y*pow((height/CurvLength),1.732);
          R[ip].x = R[ip].x;
          y[2*ip]= R[ip].x;
          y[2*ip+1] = R[ip].y;
        }
        CurvLength = 0;
        for (int ip = 1; ip < Np; ++ip){
          CurvLength = CurvLength + norm(R[ip] - R[ip-1]);
        }
      }
      // print(y,ndim);
      // exit(1);
    break;
    default:
      cout << "# We have not implemented the initial configuration:" << "\t" << niniconf << endl
           << "# EXITING" << endl;
      break;
  }
}
/****************************************************/
void iniconf_tr(double *y_tr){
  int const ndim_tr = np_tracer*pp_tracer;
  double Ymin,Ymax,Zmin,Zmax,rr,theta,radius;
  int np_tracerR,np_tracerTh,np_tracerY,left_tracer,ip_tracer;
  vec3 cc(0.02,0.0,height/2);
  switch(itracer){
    case 1:
      Ymin = -1*height/8;
      Ymax = (double)height/8;
      Zmin = 0;
      Zmax = (double)height;
      np_tracerY = (int) sqrt(np_tracer);
      for (int itracer = 0; itracer < np_tracer; ++itracer){
        y_tr[itracer*pp_tracer] = 0.01;                                         
        // X coordinate of tracer particles
      }
      for (int itracerY = 0; itracerY < np_tracerY; ++itracerY ){
        for (int itracerZ = 0; itracerZ < np_tracerY; ++itracerZ){
          ip_tracer = np_tracerY*itracerZ + itracerY;
          y_tr[ip_tracer*pp_tracer+1] =  itracerY*(Ymax-Ymin)/np_tracerY+Ymin;
          // Y coordinate of tracer particles
          y_tr[ip_tracer*pp_tracer+2] =  itracerZ*(Zmax-Zmin)/np_tracerY+Zmin;
          // Z coordinate of tracer particles
        }
      }
      left_tracer = np_tracer - np_tracerY*np_tracerY;
      for (int ip_tracer = 0; ip_tracer < left_tracer; ++ip_tracer){
        y_tr[np_tracerY*np_tracerY*pp_tracer+pp_tracer*ip_tracer+1] = (double) rand()/RAND_MAX*(Ymax-Ymin)+Ymin;
        y_tr[np_tracerY*np_tracerY*pp_tracer+pp_tracer*ip_tracer+2] = (double) rand()/RAND_MAX*(Zmax-Zmin)+Zmin;
      }
      break;
    case 2:
      rr = height/8;    // radius of the circle
      np_tracerR = (int) sqrt(np_tracer/4);   // dividing diameter and theta into equal parts.
      np_tracerTh = (int) np_tracer/np_tracerR;
      //
      for (int ip_tracer = 0; ip_tracer < np_tracer; ++ip_tracer){
        y_tr[ip_tracer*pp_tracer] = cc.x;
        // X coordinate of tracer particles
      }
      ip_tracer = 0;
      for (int itracerR = 0; itracerR < np_tracerR; ++itracerR){
        for (int itracerTh = 0; itracerTh < np_tracerTh; ++itracerTh){
          theta = 2*M_PI*itracerTh/np_tracerTh;
          y_tr[ip_tracer*pp_tracer+1] = cc.y + rr*(itracerR+1)/np_tracerR*cos(theta);
          y_tr[ip_tracer*pp_tracer+2] = cc.z + rr*(itracerR+1)/np_tracerR*sin(theta);
          ip_tracer = ip_tracer+1;
        }
      }
      left_tracer = np_tracer - np_tracerR*np_tracerTh;
      for (int itracer = 0; itracer < left_tracer; ++itracer){
        theta = 2*M_PI*rand()/RAND_MAX;
        radius = rr*rand()/RAND_MAX;
        y_tr[np_tracerR*np_tracerTh*pp_tracer+pp_tracer*itracer+1] = cc.y + radius*cos(theta);
        y_tr[np_tracerR*np_tracerTh*pp_tracer+pp_tracer*itracer+2] = cc.z + radius*sin(theta);
      }
      break;
    default:
      break;
  }
}
/********************************************/
void GetRij(vec2 R[], int i, int j, double *Distance, vec2 *rij){
  /*This calculated the distance at two index i and j on the elastic string*/
  *rij = R[j] - R[i];
  double Dis = norm(R[j]-R[i]); 
  *Distance = Dis;
}
/********************************************/
void check_param(){
  if (dd>aa){
    cout << "# ERROR: The diameter of a particle should be less than the distance between two particles." 
        << endl;
    exit(1);
  }
}
/********************************************/
void write_param( string fname ){
  ofstream paramfile;
  paramfile.open( fname );
  paramfile << "# =========== Model Parameters ==========" << endl
            << "dimension of ODE: pp = " << pp << endl
            << "Number of copies: Np = " << Np << endl
            << "height = " << height << endl
            << "aa = " << aa << endl
            << "dd = " << dd << endl
            << "viscosity = " << viscosity << endl
            << "ShearRate = " << ShearRate << endl
            << "Period = " << period << endl
            << "AA = " << AA << endl
            << "HH = " << HH << endl
            << "bcb = " << bcb << endl
            << "bct = " << bct << endl
            << "global_drag = " << UseRP << endl
            << "iext_force = " << iext_force << endl
            << "floc = " << floc << endl
            << "iext_flow = " << iext_flow << endl
            << "niniconf = " << niniconf << endl
            << "KK = " << HH*dd*dd/AA << endl
            << "MuBar = "<<  8*M_PI*viscosity*ShearRate*pow(height,4)/AA << endl;
  paramfile.close();
}
/********************************************/
void straightline(vec2 *Tngt, vec2 *Nrml, char dirn = 'y'){
  if (dirn=='y'){
    Tngt -> x = 0;
    Tngt -> y = 1;
    Nrml -> x = -1;
    Nrml -> y = 0;
  }else{
    cout<< "This function is called to set straightline in model.cpp file"
           "Only 'y' direction is implemented. EXITING!! ";
           exit(1);
  }
}
/********************************************/
vec2 y2vec2(double *y,int ip){
  vec2 XX;
  XX.x = y[2*ip];
  XX.y = y[2*ip+1];
  return XX;
}
/********************************************/
vec3 y2vec3(double *y,int ip){
  vec3 XX;
  XX.x = y[3*ip];
  XX.y = y[3*ip+1];
  XX.z = y[3*ip+2];
  return XX;
}
/********************************************/
void vec2y(double *y, vec2 XX, int ip){
  y[2*ip] = XX.x;
  y[2*ip+1] = XX.y;
}
/********************************************/
void vec2y(double *y, vec3 XX, int ip){
  y[3*ip] = XX.x;
  y[3*ip+1] = XX.y;
  y[3*ip+2] = XX.z;
}
/********************************************/
void y2kappa(double kappa[], double y[]){
  double bk, bkm1;
  // memcpy(y_start,y_start+2*pp,2*pp*sizeof(double));
  // memcpy(y_start+2*pp,y,4*sizeof(double));
  int iorbit = MM.iorbit;
  vec2 uk,ukm1;
  if(bcb==1){
    kappa[0]=0;
  }
  else{
    cout << "Boundary condition at Np==0,1 not implemented." << endl;
    cout << "Exiting" << endl;
    exit(1);
  }
  //
  if(bct==1){
    kappa[Np-1]=0;
  }else{
    cout << "Boundary condition at Np==0,1 not implemented." << endl;
    cout << "Exiting" << endl;
    exit(1);
  }
  //
  vec2 dX = y2vec2(y,1)-y2vec2(y,0);
  uk = dX/norm(dX);
  // getub(&bk,&uk,y2vec2(y,0),y2vec2(y,1));
  for (int ip = 1; ip < Np-1; ++ip){
    ukm1 = uk;
    dX = y2vec2(y,ip+1)-y2vec2(y,ip);
    uk = dX/norm(dX);
    kappa[ip] = 1/aa*cross(uk,ukm1);
  }
}
/********************************************/
void kappa2y(double y[], double kappa[]){
  // Transformation from kappa to y goes here.
  /* This is done using Frent-Serret equations.
    dT/ds = \kappa N
    dN/ds = -\kappa T
  */
  cout << "#" << y_start[0] << "\t" << y_start[1] << "\t" << y_start[2] << "\t" << y_start[3] << endl;
  double ds=aa;
  vec2 Tngt,Nrml,Tngtm1,Nrmlm1,XX,dX;
  // vec2y(y,XX,0);
  // XX.x = 0;
  // XX.y = ds;
  // vec2y(y,XX,1);
  // if (bcb==1){
  //vector from ip=0 to ip=1;
  // straightline(&Tngt,&Nrml);
  // straightline(&Tngtm1,&Nrmlm1);
  // vec2 dX = y2vec2(y_start,1)-y2vec2(y_start,0);
  // print(y_start,4);
  dX = y2vec2(y_start,1)-y2vec2(y_start,0);
  memcpy(y,y_start,2*pp*sizeof(double));
  // else if(cc=="current"){
  //   dX = y2vec2(y_start,3)-y2vec2(y_start,2);
  //   memcpy(y,y_start+2*pp,2*pp*sizeof(double));
  // }
  Tngt = dX/norm(dX);
  // Take a cross product of tngt with (\hat{k}) to get normal.
  Nrml.x = Tngt.y;
  Nrml.y = -Tngt.x;
  XX = y2vec2(y,1);
  // }
  for (int ip = 2; ip < Np; ++ip){
    Tngtm1 = Tngt;
    Nrmlm1 = Nrml;
    Tngt = Tngtm1 + Nrmlm1*kappa[ip-1]*ds;
    Tngt = Tngt/norm(Tngt);
    // Tngt = Tngt/(norm(Tngt)+tiny);
    Nrml = Nrmlm1 - Tngtm1*kappa[ip-1]*ds;
    Nrml = Nrml/norm(Nrml);
    // Nrml = Nrml/(norm(Nrml)+tiny);
    XX = XX + Tngt*ds;
    vec2y(y,XX,ip);
  }
  // Is it physical?
  bool isphysical=1;
  double d_rij;
  vec2 X1;
  vec2 X2;
  for (int ip = 0; ip < Np; ++ip){
    for (int jp = ip+1; jp < Np; ++jp){
       X1 = y2vec2(y,ip);
       X2 = y2vec2(y,jp);
       d_rij = norm(X1-X2);
       if (d_rij<dd/2){
         isphysical=0;
         break;
       }
    }
    if (isphysical==0){
      break;
    }
  }
  if (isphysical==0){
    cout << "# I am halving the curvature to take a physical NK step." << endl;
    for (int ip = 0; ip < Np; ++ip){
      kappa[ip] = kappa[ip]/2;
    }
    kappa2y(y,kappa);
  }
  // cc = "previous";
}
/********************************************/
void coordinate_transform(double y_trans[], double y[]){
  // I will write about the transformation from y to kappa here.
  y2kappa(y_trans,y);
  // memcpy(y_trans,y,ndim*sizeof(double));
  // for (int idim = 0; idim < int(ndim)/2; ++idim){
  //   y_trans[2*idim] = y_trans[2*idim] - y_trans[0];
  //   y_trans[2*idim+1] = y_trans[2*idim+1] - y_trans[1];
  // }
}
/********************************************/
void inv_coordinate_transform(double y[], double y_trans[]){
  // Transformation from kappa to y goes here.
  // string cc = "current";
  kappa2y(y,y_trans);
  // memcpy(y,y_trans,ndim*sizeof(double));
}
/********************************************/
void pre_next_iter(double *y, double *y_trans){
  cout << "# Is it coming?" << endl;
  memcpy(y_start,y,2*pp*sizeof(double));
}