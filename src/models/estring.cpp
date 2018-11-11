#include <iostream>
#include<fstream>
#include "ode.h"
#include "modules/3vec.h"
#include "modules/2Tens.h"
#include "model.h"
// #include <cmath>

/**************************/
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, double* add_SS, bool flag_kappa);
void getub(double *bk, vec3 *uk, int kp, vec3 X[]);
int MatrixtoVector(int i, int j, int N);
void GetRij(vec3 X[], int i, int j, double *Distance, vec3 *rij);
/**************************/

using namespace std;

void eval_rhs(double time,double y[],double rhs[], bool flag_kappa, double CurvSqr[], double SS[]){
  vec3 R[Np],dR[Np], EForce[Np], EForce_ip;  // R is the position of the beads.
  // double CurvSqr[Np];
  double kappasqr, MaterialCoordinate, Mobility[Np*(Np+1)/2][6];

  // Initializing Mobility Matrix.
  for (int i = 0; i < Np*(Np+1)/2; ++i)
  {
      for (int j = 0; j < 6; ++j)
      {
          Mobility[i][j] = 0;
      }
  }

  for (int ip=0;ip<Np;ip++){
    R[ip].x=y[3*ip];
    R[ip].y=y[3*ip+1];
    R[ip].z=y[3*ip+2];
  }

  for (int ip=0;ip<Np;ip++){
    kappasqr=CurvSqr[ip];
    MaterialCoordinate=SS[ip];
    dHdR(ip, R, &EForce_ip, &kappasqr, &MaterialCoordinate, flag_kappa);
    EForce[ip] = EForce_ip;
    // EForce[ip] = EForce[ip]*aa*aa;
    // dR[ip]=EForce*OneByGamma;
    CurvSqr[ip]=kappasqr;
    SS[ip]=MaterialCoordinate;
    // cout << CurvSqr[ip] << endl;
  }

  vec3 FF0(0., 0., -FFZ0*sin(omega*time));
  EForce[Np-1] = EForce[Np-1]-FF0;
  
  // cout << EForce[Np-1].x << endl;

  if (UseRP == 'O')
  {
    // cout << "Heheheheh" ;
    vec3 rij;
    double distance, c1, c2;
    for (int ip = 0; ip < Np; ++ip)
    {
      for (int jp = ip; jp < Np; ++jp)
      {
        if (jp==ip)
        {
            Mobility[MatrixtoVector(ip,ip,Np)][0] = 1.0/(3*M_PI*viscosity*dd);
            Mobility[MatrixtoVector(ip,ip,Np)][3] = Mobility[MatrixtoVector(ip,ip,Np)][0];
            Mobility[MatrixtoVector(ip,ip,Np)][5] = Mobility[MatrixtoVector(ip,ip,Np)][0];
        }
        else
        {   
            GetRij(R, ip, jp, &distance, &rij);
            c1 = 1.0/(8*M_PI*viscosity*distance)*(1+ dd*dd*1.0/(6*distance*distance));
            c2 = 1.0/(8*M_PI*viscosity*distance*distance*distance)*(1-dd*dd*1.0/(2*distance*distance));
            Mobility[MatrixtoVector(ip,jp,Np)][0] = c1+c2*(rij.x)*(rij.x);
            Mobility[MatrixtoVector(ip,jp,Np)][1] = c2*(rij.x)*(rij.y);
            Mobility[MatrixtoVector(ip,jp,Np)][2] = c2*(rij.x)*(rij.z);
            
            Mobility[MatrixtoVector(ip,jp,Np)][3] = c1+c2*(rij.y)*(rij.y);
            Mobility[MatrixtoVector(ip,jp,Np)][4] = c2*(rij.y)*(rij.z);

            Mobility[MatrixtoVector(ip,jp,Np)][5] = c1+c2*(rij.z)*(rij.z);
        }
      }
    }

    for (int ip = 0; ip < Np; ++ip)
    {
      vec3 dR_temp;
      for (int jp = 0; jp < Np; ++jp)
      {
        dR_temp.x = dR_temp.x + (Mobility[MatrixtoVector(ip,jp,Np)][0])*(EForce[jp].x)
                              + (Mobility[MatrixtoVector(ip,jp,Np)][1])*(EForce[jp].y)
                              + (Mobility[MatrixtoVector(ip,jp,Np)][2])*(EForce[jp].z);

        dR_temp.y = dR_temp.y + (Mobility[MatrixtoVector(ip,jp,Np)][1])*(EForce[jp].x)
                              + (Mobility[MatrixtoVector(ip,jp,Np)][3])*(EForce[jp].y)
                              + (Mobility[MatrixtoVector(ip,jp,Np)][4])*(EForce[jp].z);

        // if (ip == 18)
        // {
        //     cout << Mobility[MatrixtoVector(ip,jp,Np)][4] << '\t' << jp << endl;
        // }

        dR_temp.z = dR_temp.z + (Mobility[MatrixtoVector(ip,jp,Np)][1])*(EForce[jp].x)
                              + (Mobility[MatrixtoVector(ip,jp,Np)][4])*(EForce[jp].y)
                              + (Mobility[MatrixtoVector(ip,jp,Np)][5])*(EForce[jp].z);                      
      }
      dR[ip] = dR_temp;
    }

  } 

  else if (UseRP == 'Y')
  {
    // mu_ij represents the one element of mobility matrix (Size: NXN). 
    // Every element of the matrix itself is a 2nd rank tensor and the dimension of that should 3x3.

    Tens2 mu_ij, mu_ii;
    double d_rij;
    vec3 rij;

    mu_ii = dab/(3*M_PI*viscosity*dd);
    // PTens2(mu_ii);
    // rij = R[j]-R[i] and d_rij is just the norm of this value.

    for (int ip = 0; ip < Np; ++ip)
    {
        // The mu_ij in the next line represents the mobility tensor when j is equal to i and in 
        // response to that the for loop is started from ip+1 .

        dR[ip] = dR[ip] + dot(mu_ii,EForce[ip]);
        
        for (int jp = ip+1; jp < Np; ++jp)
        {
            GetRij(R, ip, jp, &d_rij, &rij);
            double c1 = 1/(8*M_PI*viscosity*d_rij);
            mu_ij = c1*(dab + rij*rij/(d_rij*d_rij) + dd*dd/(2*d_rij*d_rij)*(dab/3 - rij*rij/(d_rij*d_rij)));
            dR[ip] = dR[ip] + dot(mu_ij,EForce[jp]);
            dR[jp] = dR[jp] + dot(mu_ij,EForce[jp]);
        }
    }
  }

  else
  {
    for (int ip = 0; ip < Np; ++ip)
    {
        dR[ip] = EForce[ip]*OneByGamma;
    }
  } 
    
  
  // External force applied on the end point.
  // cout << FF0.z <<endl;
  // dR[Np-1] = dR[Np-1]-FF0*; 
  //dR[Np-1].y = 0;                     // Constraint that last point should always remain on z axis. 
  
  for (int ip=0;ip<Np;ip++){
    rhs[3*ip]=dR[ip].x;
    rhs[3*ip+1]=dR[ip].y;
    rhs[3*ip+2]=dR[ip].z;
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
void dHdR(int kp, vec3 X[], vec3* add_FF, double* add_kappasqr, double* add_SS, bool flag_kappa){
     // This function calculates the force at every node which is a function of X, time.
  vec3 ukm2(0.,0.,0.), ukm1(0.,0.,0.), uk(0.,0.,0.), ukp1(0.,0.,0.), Xzero(0.,0.,0.), dX(0.,0.,0.);
  double bkm2, bkm1, bk, bkp1;
  vec3 FF = *add_FF;             
    // Since I am passing the address of force in add_FF and the same goes for Kapppsqr

  //vec3 FF;
  Xzero.x=0.; Xzero.y=0.;Xzero.z=Z0;

  switch(kp){
    case 0:
      dX = X[kp-1+1]-Xzero;
      bkm1 = norm(dX);
      ukm1=dX/bkm1;
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      FF = (     (uk)/bkm1 - (ukm1+ukp1)/bk
           + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
           - (ukm1/bkm1)*( dot(ukm1,uk) )
           );
      FF = FF*AA/aa;
      // Add an extra term for inextensibility constraint
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa; 
      // cout << FF.z << endl;
      *add_kappasqr=0.;
      *add_FF = FF;
      *add_SS = (kp+1)*bkm1;
      break;

    case 1:
      dX = X[kp-2+1]-Xzero;
      bkm2 = norm(dX);
      ukm2 = dX/bkm2;
      getub(&bkm1, &ukm1, kp-1, X);
      getub(&bk, &uk, kp, X);
      getub(&bkp1, &ukp1, kp+1, X);
      FF = (     (uk+ukm2)/bkm1 - (ukm1+ukp1)/bk
          + (uk/bk)*( dot(uk,ukm1) + dot(uk,ukp1) )
          - (ukm1/bkm1)*( dot(ukm1,ukm2) + dot(ukm1,uk) )
          );
      FF = FF*(AA/aa);
      // cout << FF.z << endl;
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa) )*HH/aa;   // Inextensibility constraint
      *add_kappasqr=0.;
      *add_FF = FF;
      *add_SS = (kp+1)*bkm1;   
      break;

    case Np-2:
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
      *add_SS = (kp+1)*bkm1;
      break;

    case Np-1:
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
      *add_SS = (kp+1)*bkm1;
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
      // cout << bkm1 <<endl;
      // cout << FF.y << endl;
      FF = FF - (ukm1*(bkm1-aa) - uk*(bk-aa))*HH/aa;    // Inextensibility constraint 
      // cout << FF.x << endl;
      // cout << FF.y << endl;     

      if (flag_kappa)
      {
        *add_kappasqr=2.*(1.-dot(uk,ukm1));
        // cout << kappasqr << endl;
      }
      *add_FF = FF;
      *add_SS = (kp+1)*bkm1;
      break;
  }  

  //   if (flag_kappa)
  //   {
  //     *add_kappasqr=2.*(1.-dot(uk,ukm1));
  //     // cout << kappasqr << endl;
  //   }
  //   *add_FF = FF;
  // }
  // *add_SS = (kp+1)*bkm1;
 }
/**************************/
// void diagnos(int p){
  // cout << time<<'\t' ;
  // for (int ip=0;ip<Np;ip++){
  //   cout << y[3*ip] <<'\t' << y[3*ip+1]  <<'\t' << y[3*ip+2] << '\t' ;
  // }
  // cout << '\n';

  // ofstream outfile;
  // string l = "output/position";
  // l.append(to_string(itn));
  // l.append(".txt");
  // outfile.open(l, ios::out);

  // for (int ip = 0; ip < Np; ++ip)
  //   {
  //     outfile << y[3*ip] << '\t' << y[3*ip+1] << '\t' << y[3*ip+2] << endl ;
  //   }

  // outfile.close(); 
  // cout << p << endl;
  // cout << time << '\t' << y[0] << endl;
// }
/**************************/
void iniconf(double y[]){
  vec3 R[Np];  // R is the position of the beads.
  double k = 1;      // determines the frequency for initial configuration
  double CurvLength = 0;  // determines the total length of the curve
  for (int ip=0;ip<Np;ip++){
    R[ip].x=0.;
    R[ip].z=aa*double(ip+1);
    R[ip].y=aa*sin(M_PI*k*aa*double(ip+1)/height);
    // R[ip].y = 0;    
    if (ip>0)
    {
        CurvLength = CurvLength + sqrt((R[ip].x - R[ip-1].x)*(R[ip].x - R[ip-1].x)
            + (R[ip].y - R[ip-1].y)*(R[ip].y - R[ip-1].y) + (R[ip].z - R[ip-1].z)*(R[ip].z - R[ip-1].z));
        // cout << CurvLength << endl;
    }
    else
    {
        CurvLength = CurvLength + sqrt((R[ip].x)*(R[ip].x)+(R[ip].y)*(R[ip].y)+(R[ip].z)*(R[ip].z));
    }
    // cout << M_PI << endl;
    y[3*ip]=R[ip].x;
    y[3*ip+1]=R[ip].y;
    y[3*ip+2]=R[ip].z;
  }

  // R[Np-1].y = aa;
  // y[3*Np-2] = aa;
  // cout << aa << endl;
  // for (int ip = 0; ip < Np; ++ip)
  // {
  //     R[ip].y = R[ip].y/CurvLength;
  // }

}

/*********************************************/
int MatrixtoVector(int i, int j, int N)
{
  /* This function is defined to save the storage space of Mobility matrix from N^2 to N(N+1)/2 Given an 
  index of 2-D Mobility matrix this will map them to find out the same index in the 1-D Mobility vector.*/

  if (i<=j)
  {
    return (i * N - (i - 1) * i / 2 + j - i);
  }
  else
  {
    return (j * N - (j - 1) * j / 2 + i - j);  
  }

}

/********************************************/
void GetRij(vec3 R[], int i, int j, double *Distance, vec3 *rij)
{
  /*This calculated the distance at two index i and j on the elastic string*/
  // *Distance = (R[i].x - R[j].x)*(R[i].x - R[j].x) + (R[i].y - R[j].y)*(R[i].y - R[j].y) 
    // + (R[i].z - R[j].z)*(R[i].z - R[j].z);
    // The reason behind to calculate distance in this function is so that the code would waste time 
    // only once calculating the norm not every time the function norm is to be called.
    
  *rij = R[j] - R[i];
  double Dis = norm(R[j]-R[i]); 
  *Distance = Dis;
}
