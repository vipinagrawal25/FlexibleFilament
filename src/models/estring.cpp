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
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, bool flag_kappa,double time);
void getub(double *bk, vec3 *uk, int kp, vec3 X[]);
int MatrixtoVector(int i, int j, int N);
void GetRij(vec3 X[], int i, int j, double *Distance, vec3 *rij);
void drag(vec3 X[], vec3 dX[], vec3 EForce[]);
vec2 y2vec2(double *y,int ip);
vec3 y2vec3(double *y,int ip);
void ext_force(vec3* EForce,double *y,double time);
void vec2y(double *y, vec3 XX, int ip);
vec3 ext_flow(vec3 R3d, double time);
void calc_Xconstrain(vec3* Xcons, double time);
/**************************/
double y_start[2*pp] = {0.,0.,0.,0.,0.,aa};
string cc;
double Error=0;
/**************************/
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
  double EForceArr[ndim];
  eval_rhs(time,y,rhs,flag_kappa,CurvSqr,SS,EForceArr,0);
}
/**************************/
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[],
              double EForceArr, double iEforceArr){
  vec3 R[Np],dR[Np], EForce[Np], EForce_ip, FF0;  // R is the position of the beads.
  double kappasqr;
  double onebythree = 1./3.;
  SS[0] = 0;    // Initializing the material co-ordinate.
  for (int ip=0;ip<Np;ip++){
    R[ip].x=y[3*ip];
    R[ip].y=y[3*ip+1];
    R[ip].z=y[3*ip+2];
  }
  for (int ip=0;ip<Np;ip++){
    kappasqr=CurvSqr[ip];
    if (flag_kappa){
      if (ip<Np-1){
        SS[ip+1] = SS[ip] + norm(R[ip+1]-R[ip]);
      }
    }
    dHdR(ip, R, &EForce_ip, &kappasqr, flag_kappa);
    EForce[ip] = EForce_ip;
    CurvSqr[ip]=kappasqr;
  }
  if (iext_force){
    ext_force(EForce,y,time);
  }
  drag(R, dR, EForce);
  // Is the string forced at some points? Implement it here.
  for (int ip = 0; ip < Np; ++ip){
    dR[ip]=dR[ip]+ext_flow(R[ip],time);
  }
  for (int ip=0;ip<Np;ip++){
    if (iEforceArr){
      EForceArr[3*ip] = EForce[ip].x;
      EForceArr[3*ip+1] = EForce[ip].y;
      EForceArr[3*ip+2] = EForce[ip].z;
    }
    rhs[3*ip]=dR[ip].x;
    rhs[3*ip+1]=dR[ip].y;
    rhs[3*ip+2]=dR[ip].z;
  }
}
/**************************/
void ext_force(vec3* EForce, double* y, double time){
  vec3 Xghost[nconstrain],XX;
  double GG = 5000;
  if (iext_force==1 || iconstrain){
    calc_Xconstrain(Xghost,time);
    for (int i = 0; i < nconstrain; ++i){
      XX = y2vec3(y,loc_con[i]);
      EForce[i] = EForce[i] - (XX-Xghost[i])*GG;
    }
  }
}
/**************************/
vec3 ext_flow(vec3 R_ip, double time){
  double omega = 2*M_PI/period;
  vec2 dR;
  double Ts=40;
  switch(iext_flow){
    case 1:
    if (sin(omega*time)>=0){
        dR.y = dR.y + ShearRate*(height-R_ip.z)*ceil(sin(omega*time));
    }
    else{
        dR.y = dR.y + ShearRate*(height-R_ip.z)*floor(sin(omega*time));
    }
    break;
    case 2:
      // NOT SO GOOD WAY TO WRITE LIKE THIS
      if (time<Ts*(1/ShearRate)){
        dR.y = dR.y + ShearRate*(R_ip.z);
      }else{
        dR.y = dR.y + ShearRate*(R_ip.z)*sin(omega*(time-Ts*(1/ShearRate)));
      }
      break; 
    case 3:
      dR.y = dR.y + ShearRate*(height-R_ip.z)*sin(omega*time);
      break;
  }
  return dR;
}
/**************************/
void drag(vec3 X[], vec3 dX[], vec3 EForce[]){
  double onebythree = 1./3.;
  double mu0 = onebythree/(M_PI*viscosity*dd);
  if (UseRP == 'Y'){
    // mu_ij represents the one element of mobility matrix (Size: NXN). 
    // Every element of the matrix itself is a 2nd rank tensor and the dimension of that should 3x3.
    Tens2 mu_ij, mu_ii;
    double d_rij;
    vec3 rij;
    mu_ii = dab*mu0;      // dab is the unit 2 dimensional tensor. It is defined in module/2Tens file.
    for (int ip = 0; ip < Np; ++ip){
      // The mu_ij in the next line represents the mobility tensor when j is equal to i and in 
      // response to that the "for loop" is started from ip+1 .
      dX[ip] = dX[ip] + dot(mu_ii, EForce[ip]);
      for (int jp = ip+1; jp < Np; ++jp){
        GetRij(X, ip, jp, &d_rij, &rij);
        double c1 = 1/(8*M_PI*viscosity*d_rij);
        double dsqr1 = 1./(d_rij*d_rij);
        mu_ij = c1*(dab + (rij*rij)*dsqr1 + dd*dd/(2*d_rij*d_rij)*(dab*onebythree - (rij*rij)*dsqr1));
        dX[ip] = dX[ip] + dot(mu_ij, EForce[jp]);
        dX[jp] = dX[jp] + dot(mu_ij, EForce[ip]);
        // cout << dot(mu_ij, EForce[jp]).y << "\t" << dot(mu_ij, EForce[jp]).z << endl;
      }
    }
  }
  else if (UseRP == 'O'){
    Tens2 mu_ij, mu_ii;
    double d_rij;
    vec3 rij;
    mu_ii = dab*mu0;
    for (int ip = 0; ip < Np; ++ip){
      for (int jp = 0; jp < Np; ++jp){
        if (jp==ip){
          dX[ip] = dX[ip] + dot(mu_ii, EForce[ip]);
        }else{
          GetRij(X, ip, jp, &d_rij, &rij);
          double c1 = 1/(8*M_PI*viscosity*d_rij);
          double dsqr1 = 1./(d_rij*d_rij);
          mu_ij = c1*(dab + (rij*rij)*dsqr1 + dd*dd/(2*d_rij*d_rij)*(dab*onebythree - (rij*rij)*dsqr1));
          dX[ip] = dX[ip]+dot(mu_ij, EForce[jp]);
        }
      }
    }
  }
  else{
    for (int ip = 0; ip < Np; ++ip){
        dX[ip] = EForce[ip]*mu0;
    }
  } 
}
/**************************/
void calc_Xconstrain(vec3* Xcons, double time){
  switch(iconstrain){
    case 1:
      // points on circle
      for(int i = 0; i < nconstrain; ++i){
        Xcons[i].x = 0;
        Xcons[i].y = (height/2 +  aa*(loc_con[i]-loc_con[0]))*cos(angularVel*time);
        Xcons[i].z = (height/2 +  aa*(loc_con[i]-loc_con[0]))*sin(angularVel*time);
      }
      break;
    case 2:
      // points on figure '8'
      for (int i = 0; i < nconstrain; ++i){
        Xcons[i].x = 0;
        Xcons[i].y = (height +  aa*(loc_con[i]-loc_con[0]))*cos(angularVel*time);
        Xcons[i].z = (height +  aa*(loc_con[i]-loc_con[0]))*cos(angularVel*time)*sin(angularVel*time);
      }
    case 3:
      // Theta should be a triangular wave.
      double theta = 4*M_PI*abs(time/period - floor(time/period+0.5));
      for (int i = 0; i < nconstrain; ++i){
        Xcons[i].x = 0;
        Xcons[i].y = (height/2 +  aa*(loc_con[i]-loc_con[0]))*cos(theta);
        Xcons[i].z = (height/2 +  aa*(loc_con[i]-loc_con[0]))*sin(theta);
      }
    }
}
/**************************/
void calc_yone(double *yone, double time){
  double yzero[2];
  calc_yzero(yzero,time);
  yone[0] = yzero[0] + aa*cos(angularVel*time);
  yone[1] = yzero[1] + aa*sin(angularVel*time);
  // print(yone,2);
}
/**************************/
void calc_yzero(double *yzero, double time){
  yzero[0] = height/2*cos(angularVel*time);
  yzero[1] = height/2*sin(angularVel*time);
}
/**************************/
void getub(double *bk, vec3 *uk, int kp, vec3 X[]){
  vec3 dX = X[kp+1]-X[kp];
  double bb = norm(dX);
  *bk = bb;
  *uk =dX/bb;
}
/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, bool flag_kappa){
     // This function calculates the force at every node which is a function of X, time.
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.), Xzero(0.,0.,0.), dX(0.,0.,0.), Xone(0.,0.),;
  double bkm2 = 0.;
  double bkm1 = 0.;
  double bk = 0.;
  double bkp1 = 0.;
  vec2 FF = *add_FF;
  double yzero[2] = {0.,0.,0.};
  double yone[2] = {0.,0.,0.};
  vec3 FF = *add_FF;             
  // Since I am passing the address of force in add_FF and the same goes for Kapppsqr
  //vec3 FF;
  /* Here the problem is that Xzero has been taken as the first point of the rod and which is claimed to be fixed in general.
  But for some cases like the implementation of the taylor experiment, we want this to be free. For this I am implementing Xzero 
  based on the configuration. Since finally with this function we just want to calculate the force on particular node. So we can
  just change the way we calculate force on 1st and 2nd node for different configuration.*/

  /* One thing should be noted that now the expression inside case 0 and case 1 would change depending upon the configuration.
  Suppose if XZero is taken as the fixed point, then we dont need to calculate F^{0} but only F^{1} which would be implemented
  in case 0 because Xzero is the bottom most point for which we dont care to calculate the force and X[0] is the first point being
  implemented in case 0. Though for a different configuration these things should be changed.*/
  switch(kp){
    case 0:
      FF = Force_zeroth(X);
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      if (bcb==1){
        FF = ( (uk/bk)*( dot(uk,ukp1) )  - (ukp1)/bk );
        FF = FF*AA/aa;
        // Add an extra term for inextensibility constraint
        FF = FF + ( uk*(bk-aa))*HH/aa; 
        // cout << FF.z << endl;
        *add_FF = FF;
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
          dX = X[kp-2+1]-Xzero;
          bkm2 = norm(dX);
          ukm2 = dX/bkm2;
          FF = (  (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
              + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
              - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
              );
          FF = FF*(AA/aa);
          // cout << FF.z << endl;
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa) )*HH/aa;   // Inextensibility constraint
          *add_FF = FF;
      }
      else if(bcb==1){
          FF = ( (uk)/bkm1 - (ukm1+ukp1)/bk
              + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
              - (ukm1/bkm1)*( dot(ukm1,uk) )
              );
          FF = FF*(AA/aa);
          // cout << FF.z << endl;
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;   // Inextensibility constraint
          *add_FF = FF;
      }
      else{
        cout << "Boundary condition at Np==1 not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);
      }
      *add_kappasqr=0.;
      // *add_SS = (kp+1)*bkm1;      
      break;

    case Np-2:
      if(bct==1){
        getub(&bkm2, &ukm2, kp-2, X);
        getub(&bkm1, &ukm1, kp-1, X);
        getub(&bk, &uk, kp, X);
        FF = (     (uk+ukm2)/bkm1 - (ukm1)/bk
            + (uk/bk)*( dot(uk,ukm1))
            - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
            );    
        FF = FF*(AA/aa);
        // cout << bk << endl;
        // cout << FF.z << endl;
        FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
        // cout << FF.z << endl;
        *add_kappasqr=0.;
        *add_FF = FF;
      }else{
        cout << "Boundary condition at Np==NN-1 not implemented." << endl;
        cout << "Exiting" << endl;
        exit(1);

      }
      // *add_SS = (kp+1)*bkm1;
      break;

    case Np-1:
      if(bct==1){
        getub(&bkm2, &ukm2, kp-2, X);
        getub(&bkm1, &ukm1, kp-1, X);
        FF = ( (ukm2)/bkm1
          - (ukm1/bkm1)*( dot(ukm1,ukm2) )
          );
        FF = FF*(AA/aa);
        // cout << bkm1 << endl;
        // cout << FF.y << endl;
        FF = FF - (ukm1*(bkm1-aa))*HH/aa;
        // cout << FF.y << endl;
        *add_kappasqr=0.;
        *add_FF = FF;
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
      FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
      FF = FF*(AA/aa);
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
      if (flag_kappa){
        *add_kappasqr=2.*(1.- dot(uk,ukm1))/(aa*aa);
      }
      *add_FF = FF;
      break;
  }  
}
/**************************/
void iniconf(double *y){
    // Merge all three cases into 1;
  vec2 R[Np];              // R is the position of the beads.
  double k = 2;            // determines the frequency for initial configuration
  double CurvLength = 0;   // determines the total length of the curve
  string l;
  // double theta= (double)M_PI/2;
  double theta = 0;
  // double vdis=2*aa;     // Define a parameter called middle point in model.h
  double inidis=0;        // that will take care of everything
  vec2 Xcons[nconstrain];
  ifstream myfile;
  double ch;
  int cnt=0;
  switch(niniconf){
    case -1:
      myfile.open(datafile);
      // myfile >> ch;
      while(myfile >> ch){cnt++;}
      myfile.close();
      cout << "# Number of elements read: " << cnt << endl;
      if(cnt==Np){}
      else if(cnt==ndim){
        cnt=0;
        myfile.open(datafile);
        while(myfile >> ch){y[cnt]=ch;cnt++;}
        memcpy(y_start,y,2*pp*sizeof(double));
        myfile.close();
      }
      else if(cnt==pp*ndim){
        cout << datafile[0] << endl;
        myfile.open(datafile);
        for(int idim = 0; idim < Np; ++idim){
          for (int ip = 0; ip < pp; ++ip){
            myfile >> y[idim*pp+ip];
          }
          for (int ip = 0; ip < pp; ++ip){
            myfile >> ch;
          }
        }
        memcpy(y_start,y,2*pp*sizeof(double));
        myfile.close();
      }
      else{
        cout << "# Initial point is not right." << endl;
        exit(1);
      }
      myfile.close();
      break;
    case 0:
      y[0]=0; y[1]=0;y[2]=0;
      for (int ip=1;ip<Np;ip++){
        R[ip].x=0;
        R[ip].y=aa*sin(M_PI*k*aa*double(ip)/height);
        R[ip].z=aa*double(ip);
        CurvLength += norm(R[ip]-R[ip-1]);
        vec2y(y, R, ip);
      } 
      for (int idim = 0; idim < ndim; ++idim){
        y[idim] = y[idim]*height/CurvLength;
      }
      break;
    case 1:
      // In this case we implement the initial configuration for GI Taylor experiment. 
      // i.e. a straight rod which has length equal to the height of the box and free to move from bottom.
      for (int ip = 0; ip < Np; ++ip){
        R[ip].x = 0;
        R[ip].y = (aa*double(ip)+inidis)*sin(theta);
        R[ip].z = (aa*double(ip)+inidis)*cos(theta);
        //
        vec2y(y, R, ip);
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
        }
        else{
          CurvLength = CurvLength + sqrt((R[ip].x)*(R[ip].x)+(R[ip].y)*(R[ip].y));
        }
        vec2y(y,R,ip);
      }
      for (int ip = 0; ip < Np; ++ip){
        y[2*ip+1] = R[ip].y/CurvLength;
      }
      break;
    case 3:
      calc_Xconstrain(Xcons,0);
      for(int ip = 0; ip < Np; ++ip){
        R[ip].x = 0;
        R[ip].y = (aa*double(ip)+Xcons[0].x)*cos(theta);
        R[ip].z = (aa*double(ip)+Xcons[0].y)*sin(theta);
        //
        vec2y(y,R,ip);
      }
    break;
    default:
      cout << "# We have not implemented the initial configuration:" << "\t" << niniconf << endl
           << "# EXITING" << endl;
      break;
  }
}
/********************************************/
void GetRij(vec3 R[], int i, int j, double *Distance, vec3 *rij){
  /*This calculated the distance at two index i and j on the elastic string*/
  *rij = R[j] - R[i];
  double Dis = norm(R[j]-R[i]); 
  *Distance = Dis;
}
/********************************************/
void check_param(){
  // This function checks whether the parameter makes physical sense.
  if (dd>aa){
    cout << "ERROR: The diameter of a particle should be less than the distance between two particles." << endl;
    exit(1);
  }
}
/****************************************************/
void iniconf_tr(double *y_tr){
  int const ndim_tr = ntracer*ptracer;
  double Ymin = 0;
  double Ymax = (double)height/2;
  double Zmin = 0;
  double Zmax = (double)height/4;
  int ntracerY = (int) sqrt(ntracer);
  for (int itracer = 0; itracer < ntracer; ++itracer){
    y_tr[itracer*ptracer] = 0.01;                                         
    // X coordinate of tracer particles
  }
  for (int itracerY = 0; itracerY < ntracerY; ++itracerY ){
    for (int itracerZ = 0; itracerZ < ntracerY; ++itracerZ){
      int itracer = ntracerY*itracerZ + itracerY;
      y_tr[itracer*ptracer+1] =  itracerY*(Ymax-Ymin)/ntracerY+Ymin;            
      // Y coordinate of tracer particles
      y_tr[itracer*ptracer+2] =  itracerZ*(Zmax-Zmin)/ntracerY+Zmin;                   
      // Z coordinate of tracer particles
    }
  }
  int left_tracer = ntracer - ntracerY*ntracerY;
  // cout << left_tracer << endl;
  for (int itracer = 0; itracer < left_tracer; ++itracer){
    y_tr[ntracerY*ntracerY*ptracer+ptracer*itracer+1] = (double) rand()/RAND_MAX*(Ymax-Ymin)+Ymin;
    y_tr[ntracerY*ntracerY*ptracer+ptracer*itracer+2] = (double) rand()/RAND_MAX*(Zmax-Zmin)+Zmin;
  }
  // print(y_tr,ndim_tr);
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
  vec3 uk,ukm1;
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
  vec3 dX = y2vec3(y,1)-y2vec3(y,0);
  uk = dX/norm(dX);
  // getub(&bk,&uk,y2vec2(y,0),y2vec2(y,1));
  for (int ip = 1; ip < Np-1; ++ip){
    ukm1 = uk;
    dX = y2vec3(y,ip+1)-y2vec3(y,ip);
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