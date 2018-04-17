#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

int n = 101;					// total number of points on the rod.
double A = 0.001;						// A is bending rigidity constant (A = EI), unit -> Pa.m^4
double a = (double) 1/(n-1); 	// distance between two nodes, which should be constant accodrding to the constraint.
//double H =  16*A/(a);			// H should be atleast 1500 times of A. H - Stretching rigidity, unit -> Pa.m^2
double H = 0;

using namespace std;

std::vector<double> position(double p, std::vector<std::vector<double> > &R)
{
	// I have taken p as a material co-ordinate and this function returns the position of any point on the 
	// discrete rod by assuming that rod is linear in between two nearby points. It should be noted that 
	// instead of 0 to 1, p takes values from 0 to n-1.

	std::vector<double> pos(3);

	for (int i = 0; i < 3; ++i)
	{
		if (p == floor(p))
		{
			pos[i] = R[p][i];
		}
		else
		{
			pos[i] = (R[ceil(p)][i])*(p-floor(p))+(R[floor(p)][i])*(ceil(p)-p);
		}

	}

	return pos;
}

double norm(std::vector<double> &vec1, std::vector<double> &vec2)
{
	// I am here calculating the norm of a vector described as "vec1-vec2".
	double accum = 0;

	for (int i = 0; i < 3; ++i)
	{
		accum += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
	}

	accum = sqrt(accum);
	return accum;
}

double dotP(std::vector<double> &vec1, std::vector<double> &vec2, std::vector<double> &vec3)
{

	// This function computes the (vec1-vec2).(vec2-vec3);
	double accum = 0;

	for (int i = 0; i < 3; ++i)
	{
		accum += (vec1[i] - vec2[i])*(vec2[i] - vec3[i]); 
	}

	return accum;
}

double velocity(double p, int i, std::vector<std::vector<double> > &R)
{
	// This function is used to calculate velocity at every node on the road at a given time, and displacement.
	// Which is proportional to negative gradient of total energy of the system.

	std::vector<double> R2l(3), R2r(3), R1l(3), R1r(3), R0(3);
	double V, T1, T2, T3, T4, T5, T6, T7, T8;

	if (floor(p)==0)
	{
		R2r = position(p+2,R); 
		R1r = position(p+1,R);
		R0 = position(p,R);

		T4 = -(R2r[i] - R1r[i])/(norm(R2r,R1r)*norm(R1r,R0));
		T8 = (dotP(R2r,R1r,R0)*(R1r[i]-R0[i]))/(norm(R2r,R1r)*pow(norm(R1r,R0),3.0));	

		// V = -A*(R2r[i]-R1r[i])/(a*a*a) 
		V =  A/a*(T4+T8) + H*(a*(R0[i]-R1r[i])/norm(R1r,R0) + R1r[i]- R0[i])/a;	
		// cout << H*(a*(R0[i]-R1r[i])/norm(R1r,R0) + R1r[i]- R0[i])/a << endl;
	}
	else if (floor(p) == 1)
	{
		R2r = position(p+2,R);
		R1r = position(p+1,R);
		R1l = position(p-1,R);
		R0 = position(p,R);

		T2 = -(R0[i] - R1l[i])/(norm(R1r,R0)*norm(R0,R1l));
		T3 = (R1r[i] - R0[i])/(norm(R1r,R0)*norm(R0,R1l));
		T4 = -(R2r[i] - R1r[i])/(norm(R2r,R1r)*norm(R1r,R0));
		T6 = (dotP(R1r,R0,R1l)*(R1r[i]-R0[i]))/(norm(R0,R1l)*pow(norm(R1r,R0),3.0));		
		T7 = -(dotP(R1r,R0,R1l)*(R0[i]-R1l[i]))/(norm(R1r,R0)*pow(norm(R0,R1l),3.0));		
		T8 = (dotP(R2r,R1r,R0)*(R1r[i]-R0[i]))/(norm(R2r,R1r)*pow(norm(R1r,R0),3.0));

		// V = -A*(R2r[i]-2*R1r[i]+2*R0[i]-R1l[i])/(a*a*a) 
		// cout << "I am here" << norm(R1r,R0) << endl;
		V = A/a*(T2+T3+T4+T6+T7+T8)- H*(2*R0[i]-R1l[i]-R1r[i]-a*(R0[i] - R1l[i])/norm(R0,R1l) - a*(R0[i] - R1r[i])/norm(R0,R1r))/a;
		// cout << H*(2*R0[i]-R1l[i]-R1r[i]-a*(R0[i] - R1l[i])/norm(R0,R1l) - a*(R0[i] - R1r[i])/norm(R0,R1r))/a << endl;
	}
	else if (ceil(p) == n-2)
	{
		R1r = position(p+1,R);
		R1l = position(p-1,R);
		R2l = position(p-2,R);
		R0 = position(p,R);

		T1 = (R1l[i] - R2l[i])/(norm(R0,R1l)*norm(R1l,R2l));
		T2 = -(R0[i] - R1l[i])/(norm(R1r,R0)*norm(R0,R1l));
		T3 = (R1r[i] - R0[i])/(norm(R1r,R0)*norm(R0,R1l));
		T5 = -(dotP(R0,R1l,R2l)*(R0[i]-R1l[i]))/(norm(R1l,R2l)*pow(norm(R0,R1l),3.0));
		T6 = (dotP(R1r,R0,R1l)*(R1r[i]-R0[i]))/(norm(R0,R1l)*pow(norm(R1r,R0),3.0));
		T7 = -(dotP(R1r,R0,R1l)*(R0[i]-R1l[i]))/(norm(R1r,R0)*pow(norm(R0,R1l),3.0));

		// cout << "This is useful too" << endl;
		// V = A*(R1r[i]+2*R1l[i]-R2l[i]-2*R0[i])/(a*a*a)
		V = A/a*(T1+T2+T3+T5+T6+T7)- H*(2*R0[i]-R1l[i]-R1r[i]-a*(R0[i] - R1l[i])/norm(R0,R1l) - a*(R0[i] - R1r[i])/norm(R0,R1r))/a;
		//cout << i << '\t' << H*(2*R0[i]-R1l[i]-R1r[i]-a*(R0[i] - R1l[i])/norm(R0,R1l) - a*(R0[i] - R1r[i])/norm(R0,R1r))/a << endl;
		cout << T1 << '\t' << T2 << '\t' << T3 << '\t' << T5 << '\t' << T6 << '\t' << T7 << endl;
	}
	else if (ceil(p)==n-1)
	{
		R2l = position(p-2,R); 
		R1l = position(p-1,R);
		R0 = position(p,R);

		T1 = (R1l[i] - R2l[i])/(norm(R0,R1l)*norm(R1l,R2l));
		T5 = -(dotP(R0,R1l,R2l)*(R0[i]-R1l[i]))/(norm(R1l,R2l)*pow(norm(R0,R1l),3.0));

		//cout << "Why the fuck loop is not coming here?" << endl;
		// V = A*(R1l[i]-R2l[i])/(a*a*a) 
		V = A/a*(T1+T5)+ H*(a*(R0[i]-R1l[i])/norm(R0,R1l) +R1l[i]- R0[i])/a;
		//cout << V << endl;
	}
	else		
	{	
		R2l = position(p-2,R);
		R2r = position(p+2,R);
		R1l = position(p-1,R);
		R1r = position(p+1,R);
		R0 = position(p,R);


		T1 = (R1l[i] - R2l[i])/(norm(R0,R1l)*norm(R1l,R2l));
		T2 = -(R0[i] - R1l[i])/(norm(R1r,R0)*norm(R0,R1l));
		T3 = (R1r[i] - R0[i])/(norm(R1r,R0)*norm(R0,R1l));
		T4 = -(R2r[i] - R1r[i])/(norm(R2r,R1r)*norm(R1r,R0));
		T5 = -(dotP(R0,R1l,R2l)*(R0[i]-R1l[i]))/(norm(R1l,R2l)*pow(norm(R0,R1l),3.0));
		T6 = (dotP(R1r,R0,R1l)*(R1r[i]-R0[i]))/(norm(R0,R1l)*pow(norm(R1r,R0),3.0));
		T7 = -(dotP(R1r,R0,R1l)*(R0[i]-R1l[i]))/(norm(R1r,R0)*pow(norm(R0,R1l),3.0));
		T8 = (dotP(R2r,R1r,R0)*(R1r[i]-R0[i]))/(norm(R2r,R1r)*pow(norm(R1r,R0),3.0));

		// V = A*(-R2l[i] + 2*R1l[i] - 2*R0[i] + 2*R1r[i] -R2r[i])/(a*a*a)
		V = A/a*(T1+T2+T3+T4+T5+T6+T7+T8)- H*(2*R0[i]/a-R1l[i]/a-R1r[i]/a-(R0[i] - R1l[i])/norm(R0,R1l) - (R0[i] - R1r[i])/norm(R0,R1r));
		// cout << T1 << '\t' << T2 <<'\t' << T3<< endl;
	}

	return V;	
}

int main()
{
	// double V0 = 0.1;	// V0 is the velocity of two end points of the rod in z direction.	
	int steps = 1;
	double dt = 0.01;
	double p = 0;
	double F = 0.001;		// F is the force which is acting on the rod.								

	std::vector<std::vector<double> > R(n, std::vector<double> (3));	// R is the position of the beads.
	std::vector<std::vector<double> > R_temp(n, std::vector<double> (3));

	double k1,k2,k3,k4;

	for (int i = 0; i < n; ++i)
	{
		R[i][0] = 0;
		R[i][1] = 0;
		R[i][2] = (double) i/(n-1);
	}

	R_temp = R;
	while (steps<800)
	{	
		// R[0][2] = R[0][2] + V0*dt;
		// R[n-1][2] = R[n-1][2] - V0*dt;

		for (int i = 1; i < n-1; ++i)
		{
			p = (double) i;

			for (int j = 0; j < 3; ++j)
			{
				k1 = velocity(p,j,R);
				//cout << floor(p) << '\t' << k1 << endl;				
				k2 = velocity(p+dt*k1/2,j,R);	
				k3 = velocity(p+dt*k2/2,j,R);
				k4 = velocity(p+dt*k3,j,R);
				R_temp[i][j] = R[i][j] + dt*(k1+2*k2+2*k3+k4)/6;
				//cout << floor(p) << '\t' << dt*(k1+2*k2+2*k3+k4)/6 << endl;
			}
		}

		R_temp[0][2] = 0;
		R_temp[n-1][2] = R_temp[n-1][2] - F*dt;
		R = R_temp;
		// cout << k1 << '\t' << k2 << '\t' << k3 << endl;
		ofstream outfile;
	    outfile.open("position.txt", ios::out); 
	    
		for (int i = 0; i < n; ++i)
		{
			outfile << R[i][0] << ' ' << R[i][1] << ' ' << R[i][2] << endl;
		}
	
	    outfile.close();

	    steps++;

	    cout << steps << endl;
	}
	return 0;
}
