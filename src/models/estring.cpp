#include <iostream>
#include <fstream>
#include "ode.h"
#include "modules/3vec.h"
#include "modules/2Tens.h"
#include "model.h"
#include<string>
#include<vector>
#include <math.h>
using namespace std;
/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, bool flag_kappa);
void getub(double *bk, vec3 *uk, int kp, vec3 X[]);
int MatrixtoVector(int i, int j, int N);
void GetRij(vec3 X[], int i, int j, double *Distance, vec3 *rij);
void drag(vec3 X[], vec3 dX[], vec3 EForce[]);
/**************************/
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
vec3 R[Np],dR[Np], EForce[Np], EForce_ip, FF0;  // R is the position of the beads.
// double CurvSqr[Np];
double kappasqr;
double onebythree = 1./3.;
SS[0] = 0;    // Initializing the material co-ordinate
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
drag(R, dR, EForce);
// Is the string forced at some points? Implement it here.
if (iext_force){
  FF0.x = 0;
  FF0.y = 0;
  FF0.z = -FFZ0*sin(omega*time);
  if (floc>=0 && floc<Np){
    EForce[floc] = EForce[Np-1]-FF0;
  }
  else{
    cout << "Location of the external force is wrong." <<endl;
    cout << "Exiting" << endl;
    exit(1); 
  }
}
switch(iext_flow){
   case 1:
    if (sin(omega*time)>=0){
      for (int ip = 0; ip < Np; ++ip){
        dR[ip].y = dR[ip].y + ShearRate*(height-R[ip].z)*ceil(sin(omega*time));
      }
    }
    else{
       for (int ip = 0; ip < Np; ++ip){
        dR[ip].y = dR[ip].y + ShearRate*(height-R[ip].z)*floor(sin(omega*time));
      }
    }
    break;
  case 2:
    for (int ip = 0; ip < Np; ++ip){
      double Ts = 40;
      // NOT SO GOOD WAY TO WRITE LIKE THIS
      if (time<Ts*(1/ShearRate)){
        dR[ip].y = dR[ip].y + ShearRate*(R[ip].z);
      }else{
        dR[ip].y = dR[ip].y + ShearRate*(R[ip].z)*sin(omega*(time-Ts*(1/ShearRate)));
      }
    }
    break; 
  case 3:
    for(int ip=0; ip<Np;ip++){
      dR[ip].y = dR[ip].y + ShearRate*(height-R[ip].z)*sin(omega*time);
    }
    break;
}
for (int ip=0;ip<Np;ip++){
  rhs[3*ip]=dR[ip].x;
  rhs[3*ip+1]=dR[ip].y;
  rhs[3*ip+2]=dR[ip].z;
  }
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
    mu_ii = dab*mu0;    // dab is the unit 2 dimensional tensor. It is defined in module/2Tens file.
    // PTens2(mu_ii);
    // rij = R[j]-R[i] and d_rij is just the norm of this value.

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
void getub(double *bk, vec3 *uk, int kp, vec3 X[]){
  vec3 dX = X[kp+1]-X[kp];
  double bb = norm(dX);
  *bk = bb;
  *uk =dX/bb;
}
/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, bool flag_kappa){
     // This function calculates the force at every node which is a function of X, time.
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.), Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec3 FF = *add_FF;             
  // Since I am passing the address of force in add_FF and the same goes for Kapppsqr
  //vec3 FF;
  
  if (bcb==0){
    // If the bottom point is fixed, then it can not move.
    Xzero.x=0.; Xzero.y=0.; Xzero.z=Z0;      
  }
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
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      if (bcb==0){
          dX = X[kp-1+1]-Xzero;
          bkm1 = norm(dX);
          // cout << bkm1 << endl;
          ukm1=dX/bkm1;
          FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
               + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
               - (ukm1/bkm1)*( dot(ukm1,uk) )
               );
          FF = FF*AA/aa;
          // Add an extra term for inextensibility constraint
          FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa; 
          // cout << FF.z << endl;
          *add_FF = FF;
          // *add_SS = (kp+1)*bkm1;
      }
      else if(bcb==1){
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
    
        FF = (     (ukm2)/bkm1
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
      if (flag_kappa==false){
        *add_kappasqr=2.*(1.- dot(uk,ukm1))/(aa*aa);
      }
      *add_FF = FF;
      break;
  }  
}
/**************************/
void iniconf(double *y, double *vel){
    vec3 R[Np];              // R is the position of the beads.
    double k = 1;            // determines the frequency for initial configuration
    double CurvLength = 0;   // determines the total length of the curve
    if (lastfile){
      string l = "output/var";
      l.append(to_string(lastfile));
      l.append(".txt");
      ifstream myfile(l,ios::in); 
      double num = 0.0;           
      for (int ip = 0; ip < Np; ++ip){
        myfile >> y[3*ip];
        myfile >> y[3*ip+1];
        myfile >> y[3*ip+2];
        // Now just throw next three numbers as they contain values of velocity.
        myfile >> vel[3*ip];
        myfile >> vel[3*ip+1];
        myfile >> vel[3*ip+2];
      }
      myfile.close();  
    }
    else{
      switch(niniconf){
        case 0:
          for (int ip=0;ip<Np;ip++){
            R[ip].x=0.;
            R[ip].y=aa*sin(M_PI*k*aa*double(ip+1)/height);
            R[ip].z=aa*double(ip+1);  
            // R[ip].y = 0;    
            if (ip>0){
                CurvLength = CurvLength + norm(R[ip]-R[ip-1]);
                // cout << CurvLength << endl;
            }else{
                CurvLength = CurvLength + sqrt((R[ip].x)*(R[ip].x)+(R[ip].y)*(R[ip].y)+(R[ip].z)*(R[ip].z));
            }
            y[3*ip]=R[ip].x;
            y[3*ip+1]=R[ip].y;
            y[3*ip+2]=R[ip].z;

            // Compute velocity here or in solve.cpp. It is just an evaluation of eval_rhs function.
          }
          break;
        case 1:
          // In this case we implement the initial configuration for GI Taylor experiment. 
          // i.e. a straight rod which has length equal to the height of the box and free to move from bottom.
          for (int ip = 0; ip < Np; ++ip){
              R[ip].x = 0;
              R[ip].y = 0;
              R[ip].z = aa*double(ip);
              // cout << R[ip].z << endl ;
              y[3*ip] = R[ip].x;
              y[3*ip+1] = R[ip].y;
              y[3*ip+2] = R[ip].z;
          }
          break;
        case 2:
          // Santillian's experiment
          // In this case, we want to study the dynamics of a rod which is kept in the direction of the flow at origin. The rod
          // should be deviated a little bit from origin in starting.           
          for (int ip = 0; ip < Np; ++ip)
          {
              R[ip].x = 0;
              R[ip].y = aa*(double(ip+1))-height*0.5;
              R[ip].z = aa*sin(M_PI*k*aa*double(ip+1)/height);
              if (ip>0){
                CurvLength = CurvLength + norm(R[ip] - R[ip-1]);
                  // cout << CurvLength << endl;
              }
              else{
                CurvLength = CurvLength + sqrt((R[ip].x)*(R[ip].x)+(R[ip].y)*(R[ip].y)+(R[ip].z)*(R[ip].z));
              }
              y[3*ip] = R[ip].x;
              y[3*ip+1] = R[ip].y;
              y[3*ip+2] = R[ip].z;
          }
          for (int ip = 0; ip < Np; ++ip){
              R[ip].z= R[ip].z/CurvLength;
              y[3*ip+2] = R[ip].z;
          }
          break;
      }
    }
}

/********************************************/
void GetRij(vec3 R[], int i, int j, double *Distance, vec3 *rij)
{
  /*This calculated the distance at two index i and j on the elastic string*/

  *rij = R[j] - R[i];
  double Dis = norm(R[j]-R[i]); 
  *Distance = Dis;
}
/********************************************/
void check_param(double dd, double aa){
if (dd>aa){
      cout << "ERROR: The diameter of a particle should be less than the distance between two particles." << endl;
      exit(0);
  }
}