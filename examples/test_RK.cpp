#include <iostream>

#include "RungeKutta.hpp"

using namespace std;


// ---------------------------------------------------------------------------- //

int main()
{
  double Z = 8.;
  double lambda = 1.;
  
  RungeKutta RK(1,1,6);

/*
  double sol = exp(lambda*Z);
  double err[10];
  
  for(unsigned int j = 0 ; j < 10 ; j++) {
    unsigned int NZ = 8 * pow(2,j);
    double h = Z/NZ;
  
    double y = 1.;
    double y_next;
  
    for(unsigned int i = 0 ; i < NZ ; i++) {
      y_next = RK.apply_method(h,lambda,y);
      y = y_next;
    }
    
    err[j] = abs(y-sol);
  }
  */

  double sol = exp(exp(lambda*sin(Z))-1.);

  double err[10];
  for(unsigned int j = 0 ; j < 10 ; j++) {
    unsigned int NZ = 8 * pow(2,j);
    double h = Z/NZ;
  
    double y = 1.;
    double y_next;
    
    double z = 0.;

    for(unsigned int i = 0 ; i < NZ ; i++) {
      y_next = RK.apply_method(h,z,lambda,y);
      y = y_next;
      z += h;
    }

    err[j] = abs(y-sol);
  }
  
  cout<<"err = [";
  for(unsigned int j = 0 ; j < 10 ; j++) {
    cout<<err[j]<<",";
  }
  cout<<"];"<<endl;
  
  cout<<"delta_z = [";
  for(unsigned int j = 0 ; j < 10 ; j++) {
    unsigned int NZ = 8 * pow(2,j);
    double h = Z/NZ;
    cout<<h<<",";
  }
  cout<<"];"<<endl;
  
  return 0;
}
