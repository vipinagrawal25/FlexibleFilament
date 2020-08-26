#include <iostream>
#include <fstream>
#include "ode.h"
#include "modules/2vec.h"
#include "modules/2b2Tens.h"
#include "model.h"
#include <string>
#include <math.h>
#include <sstream>
#include "constant.h"
#include <sys/stat.h>
/**************************/
using namespace std;
/**************************/
/* 1) We consider simulation of elastic filament in 2 dimension. Filament is in XY dimension.
      X direction is horizontal. */
/**************************/
void dHdR(int kp, vec2 X[], vec2* add_FF, double* add_kappasqr, bool flag_kappa);
void getub(double *bk, vec2 *uk, int kp, vec2 X[]);
int MatrixtoVector(int i, int j, int N);
void GetRij(vec2 X[], int i, int j, double *Distance, vec2 *rij);
void drag(vec2 X[], vec2 dX[], vec2 EForce[]);
bool IsPathExist(const std::string &s);
/**************************/
void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
vec2 R[Np],dR[Np],EForce[Np],EForce_ip,FF0;  
// R is the position of the beads.
// double CurvSqr[Np];
double kappasqr;
double onebythree = 1./3.;
SS[0] = 0;    // Initializing the material co-ordinate
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
dHdR(ip, R, &EForce_ip, &kappasqr, flag_kappa);
EForce[ip] = EForce_ip;
CurvSqr[ip]=kappasqr;
}
drag(R, dR, EForce);
// Is the string forced at some points? Implement it here.
if (iext_force){
  FF0.x = 0;
  FF0.y = -FFY0*sin(omega*time);
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
        dR[ip].x = dR[ip].x + ShearRate*(height-R[ip].y)*ceil(sin(omega*time));
      }
    }
    else{
       for (int ip = 0; ip < Np; ++ip){
        dR[ip].x = dR[ip].x + ShearRate*(height-R[ip].y)*floor(sin(omega*time));
      }
    }
    break;
  case 2:
    for (int ip = 0; ip < Np; ++ip){
      double Ts = 40;
      // NOT SO GOOD WAY TO WRITE LIKE THIS
      if (time<Ts*(1/ShearRate)){
        dR[ip].x = dR[ip].x + ShearRate*(R[ip].y);
      }else{
        dR[ip].x = dR[ip].x + ShearRate*(R[ip].y)*sin(omega*(time-Ts*(1/ShearRate)));
      }
    }
    break; 
  case 3:
    for(int ip=0; ip<Np;ip++){
      dR[ip].x = dR[ip].x + ShearRate*(height-R[ip].y)*sin(omega*time);
    }
    break;
}
for (int ip=0;ip<Np;ip++){
  rhs[2*ip]=dR[ip].x;
  rhs[2*ip+1]=dR[ip].y;
  }
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
    mu_ii = dab2b2*mu0;    // dab2b2 is the unit 2 dimensional tensor. It is defined in module/2b2Tens file.
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
            mu_ij = c1*(dab2b2 + (rij*rij)*dsqr1 + dd*dd/(2*d_rij*d_rij)*(dab2b2*onebythree - (rij*rij)*dsqr1));
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
void dHdR(int kp, vec2 X[], vec2* add_FF, double* add_kappasqr, bool flag_kappa){
  // This function calculates the force at every node which is a function of X, time.
  vec2 ukm2(0.,0.), ukm1(0.,0.), uk(0.,0.), ukp1(0.,0.), Xzero(0.,0.), dX(0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec2 FF = *add_FF;             
  // Since I am passing the address of force in add_FF and the same goes for Kapppsqr
  //vec2 FF;
  if (bcb==0){
    // If the bottom point is fixed, then it can not move.
    Xzero.x=0.; Xzero.y=0.;      
  }
  /* Here the problem is that Xzero has been taken as the first point of the rod and which is claimed to be fixed in general.
  But for some cases like the implementation of the taylor experiment, we want this to be free. For this I am implementing 
  Xzero based on the configuration. Since finally with this function we just want to calculate the force on particular node. 
  So we can just change the way we calculate force on 1st and 2nd node for different configuration.*/

  /* One thing should be noted that now the expression inside case 0 and case 1 would change depending upon the configuration.
  Suppose if XZero is taken as the fixed point, then we dont need to calculate F^{0} but only F^{1} which would be implemented
  in case 0 because Xzero is the bottom most point for which we dont care to calculate the force and X[0] is the first point 
  being implemented in case 0. Though for a different configuration these things should be changed.*/
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
/****************************************************/
void iniconf(double *y){
    vec2 R[Np];              // R is the position of the beads.
    double k = 1;            // determines the frequency for initial configuration
    double CurvLength = 0;   // determines the total length of the curve
    string l;
    ifstream myfile;
    double theta=ini_theta;
    if (lastfile){
      switch(wDataMeth){
        case 1:
          if (IsPathExist("output")){l = "output/var";}else{l="var";}
          l.append(to_string(lastfile));
          l.append(".txt");
          myfile.open(l,ifstream::in); 
          rData(&myfile,&y[0]);
          myfile.close();
          break;
        case 2:
          if (IsPathExist("data")){l = "data/PSI";}else{l="PSI";}
          // l.append(to_string(lastfile));
          myfile.open(l,ifstream::in);
          rData(&myfile,&y[0]);
          myfile.close();
          break;
        default:
          cout << "Reading method not implemented." << endl;
          exit(1);
      }
    }
    else{
      switch(niniconf){
        case 0:
          for (int ip=0;ip<Np;ip++){
            R[ip].x=aa*sin(M_PI*k*aa*double(ip+1)/height);
            R[ip].y=aa*double(ip+1);
            // R[ip].y = 0;    
            if (ip>0){
                CurvLength = CurvLength + norm(R[ip]-R[ip-1]);
                // cout << CurvLength << endl;
            }else{
                CurvLength = CurvLength + sqrt((R[ip].x)*(R[ip].x)+(R[ip].y)*(R[ip].y));
            }
            y[2*ip]=R[ip].x;
            y[2*ip+1]=R[ip].y;
            // Compute velocity here or in solve.cpp. It is just an evaluation of eval_rhs function.
            for (int ip = 0; ip < Np; ++ip){
              y[2*ip+1] = R[ip].y/CurvLength;
            }
          }
          break;
        case 1:
          // In this case we implement the initial configuration for GI Taylor experiment. 
          // i.e. a straight rod which has length equal to the height of the box and free to move from bottom.
          for (int ip = 0; ip < Np; ++ip){
              R[ip].x = aa*double(ip)*sin(theta);
              R[ip].y = aa*double(ip)*cos(theta);
              // cout << R[ip].z << endl ;
              y[2*ip] = R[ip].x;
              y[2*ip+1] = R[ip].y;
          }
          break;
        case 2:
          // Santillian's experiment
          // In this case, we want to study the dynamics of a rod which is kept in the direction of the flow at origin. 
          // The rod should be deviated a little bit from origin in starting.           
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
        // case 3:
        //   // The filament is rotated by some angle. theta=Pi/2 would mean a straight filament
        //   // in the direction of the flow.
        //   for (int ip = 0; ip < Np; ++ip){
        //       R[ip].x = aa*double(ip)*cos(theta);
        //       R[ip].y = aa*double(ip)*sin(theta);
        //       // cout << R[ip].z << endl ;
        //       y[2*ip] = R[ip].x;
        //       y[2*ip+1] = R[ip].y;
        //   }
          break;
      }
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
void rData(ifstream *fptr, double *y){
  string line,token;
  double num;
  istringstream iss;
  switch(rDataMeth){
    case 1:
      num=0.0;
      for(int ip = 0; ip < Np; ++ip){
        *fptr >> y[2*ip];
        *fptr >> y[2*ip+1];
        // Now just throw away next two numbers as they contain values of velocity.
        *fptr >> num;
        *fptr >> num;
      }
      break;
    case 2:
      while( getline(*fptr,token) ){
        line=token;
        getline(*fptr,token);          // Dumping this: #---------
      }
      // Now convert all the tab separated entries to array.
      iss.str(line);
      getline(iss, token, '\t');        // First entry is time, delete that.
      for (int idim = 0; idim < ndim; ++idim){
        getline(iss, token, '\t');       // Get next token.
        y[idim]=stod(token);            // Convert to double and store it in y.
      }
      break;
    default:
      cout << "Hey, your choice of writing data does not exist. "
              "If you want a new way, Code function: wData(ofstream *fptr, double y[], double vel[]) "
              "in model.cpp file." << endl;
      exit(1);
  }
}
/********************************************/
void check_param(){
  if (dd>aa){
    cout << "ERROR: The diameter of a particle should be less than the distance between two particles." << endl;
    exit(1);
  }
}
/********************************************/
void write_param( string fname ){
  ofstream paramfile;
  paramfile.open( fname );
  paramfile << "# =========== Model Parameters ==========" << endl
            << "#dimension of ODE: pp = " << pp << endl
            << "#Number of copies: Np = " << Np << endl
            << "height = " << height << endl
            << "aa = " << aa << endl
            << "dd = " << dd << endl
            << "viscosity = " << viscosity << endl
            << "Z0 = " << Z0 << endl
            << "Famp = " << FFY0 << endl
            << "sigma = " << sigma << endl
            << "ShearRate = " << ShearRate << endl
            << "omega = " << omega << endl
            << "AA = " << AA << endl
            << "HH = " << HH << endl
            << "bcb = " << bcb << endl
            << "bct = " << bct << endl
            << "global_drag = " << UseRP << endl
            << "iext_force = " << iext_force << endl
            << "floc = " << floc << endl
            << "iext_flow = " << iext_flow << endl
            << "iniconf = " << iniconf << endl
            << "KK = " << HH*dd*dd/AA << endl
            << "Gamma = "<<  8*M_PI*viscosity*ShearRate*dd*dd*dd*height/AA << endl;
  paramfile.close();  
}
/********************************************/
// Move it to utilities
bool IsPathExist(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}
/********************************************/
void reduceSymmetry(double y[]){
  // 1st is horizontal continuous symmetry
  // 2nd is vertical continuous symmetry 
  // 3rd is mix of previous two
  // 4th is a discrete symmetry if the system is rotated by an angle of Pi/2
}