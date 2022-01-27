#include <iostream>

#include "LawsonRK.hpp"

using namespace std;

typedef vector<ComplexArray> MultipleComplexArrays;
typedef vector<double> doubleArray;

// ---------------------------------------------------------------------------- //

int main()
{
  unsigned int M = 1;

  unsigned int NT = (int) pow(2,10);
  double T = 2*M_PI;

  vector<doubleArray> beta(3,doubleArray(3,1.));

  vector<Sparse3DMatrix> Q(1);
  Q[0].addNonZeroValue(0,0,0,1.);

  double Z = 2*M_PI;
  unsigned int NZ = 64;

  unsigned int order_RK = 6;

  double gamma = 1.;

  double fR = 0.;
  double tau[2] = {0.,0.};

  MultipleComplexArrays phi_in(1,ComplexArray(NT));
  
  complex<double> i(0.,0.);
  double kp = 2*M_PI;
  double t = 0.;
  double delta_t = T/NT;
  for(unsigned int j = 0 ; j < NT ; j++) {
    phi_in[0][j] = exp(i*(kp*t/T));
    // phi_0[j] = complex<double>(1./(2+cos(kp/T*t)),0.);

    t += delta_t;
  }
  
  MultipleComplexArrays phi_out(1,ComplexArray(NT));
  
  double delta_z = Z/NZ;
  LawsonRK solver(M,1.,NT,T,order_RK,1,beta,Q,fR,tau);
  
  solver.initializeLawson(phi_in,delta_z);
  solver.compute(delta_z,NZ,Z);
  
  ComplexArray phi_petit(NT);
  phi_petit = phi_out[0];

  NZ *= 2;
  delta_z /= 2;
  solver.initializeLawson(phi_in,delta_z);
  solver.compute(delta_z,NZ,Z);
  ComplexArray phi(NT);
  phi = phi_out[0];
  
  NZ *= 2;
  delta_z /= 2;
  solver.initializeLawson(phi_in,delta_z);
  solver.compute(delta_z,NZ,Z);
  ComplexArray phi_grand(NT);
  phi_grand = phi_out[0];
  
  ComplexArray num(NT);
  num = phi - phi_grand;
  ComplexArray den(NT);
  den = phi_petit - phi;
  
  cout<<"norm(num) = "<<norm(num)<<endl;
  cout<<"norm(den) = "<<norm(den)<<endl;
  
  double p = -log2(norm(num)/norm(den));
  cout<<"p = "<<p<<endl;
  
  return 0;
}
