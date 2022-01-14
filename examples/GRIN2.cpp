#include <iostream>

#include "MultimodePropagation.hpp"
#include "writeMatFile.hpp"

using namespace std;

typedef vector<ComplexArray> MultipleComplexArray;
typedef vector<double> doubleArray;

ComplexArray build_soliton(double t_final, unsigned int nt, double t0, double c, double q);

// ---------------------------------------------------------------------------- //

int main()
{
  unsigned int nt = (int) pow(2,20);
  double t_final = 160.; // in ps

  double t0 = -t_final/10;
  double c = 1.;

  // ----------- Set up input parameters ----------- //
  
  struct SimulationParameters in;
  
  in.n_modes = 2;

  vector<double> beta_0({0., 0.});
  vector<double> beta_1({0., 0.});
  vector<double> beta_2({18.938e-3, 18.866e-3});

  in.dispersion_coefficients.push_back(beta_0);
  in.dispersion_coefficients.push_back(beta_1);
  in.dispersion_coefficients.push_back(beta_2);
  
  double SR = 6.4341e12;
  
  Sparse3DMatrix Q0;  
  Q0.addNonZeroValue(0,0,0,SR);
  // Q0.addNonZeroValue(0,0,1,SR/3);
  // Q0.addNonZeroValue(0,1,1,SR/2);
  // // Q0.addNonZeroValue(1,0,0,SR/2);
  // // Q0.addNonZeroValue(1,1,0,SR/3);
  // Q0.addNonZeroValue(1,1,1,-SR/3);
  in.coupling_coefficients.push_back(Q0);
  
  Sparse3DMatrix Q1;
  // Q1.addNonZeroValue(0,0,0,-SR/3);
  // // Q1.addNonZeroValue(0,0,1,SR/3);
  // // Q1.addNonZeroValue(0,1,1,SR/2);
  // Q1.addNonZeroValue(1,0,0,SR/5);
  // Q1.addNonZeroValue(1,1,0,SR/3);
  Q1.addNonZeroValue(1,1,1,0.8*SR);
  in.coupling_coefficients.push_back(Q1);
  
  // in.fiber_length = t_final / ((sqrt(beta_2[0]/2))) / 4; // in m
  in.fiber_length = 32.;
  cout<<"Fiber length = "<<in.fiber_length<<endl;

  in.n_steps = (int) (65536/2);

  in.order_RK = 6;

  double lambda = 1030e-9; // central wavelength in m
  double f0 = 2.99792458e-4/lambda;
  double speed_of_light = 2.99792458e-4; // speed of light in m/ps
  double w0 = 2*M_PI*f0;
  double n2 = 2.3*1e-20; // m^2 W^-1
  in.nonlinearity_const = n2*w0/speed_of_light;

  in.raman_proportion = 0.;
  in.raman_parameters[0] = 12.2e-3; // in ps
  in.raman_parameters[1] = 32.e-3; // in ps

  in.nt = nt;
  in.pulse_width = sqrt(beta_2[0])/(SR*in.nonlinearity_const);
  in.time_window = t_final;

  MultipleComplexArray phi_0(2,ComplexArray(nt));
  phi_0[0] = build_soliton(t_final/sqrt(beta_2[0]/2),nt,t0/sqrt(beta_2[0]/2),c*sqrt(2/beta_2[0]),SR*in.nonlinearity_const);
  phi_0[1] = build_soliton(t_final/sqrt(beta_2[0]/2),nt,-t0/sqrt(beta_2[0]/2),-c*sqrt(2/beta_2[0]),0.8*SR*in.nonlinearity_const);
  
  for(unsigned int p = 0 ; p < 2 ; p++) {
    in.signal.push_back(phi_0[p]);
  }

  MultimodePropagation mm_propagation(in);
  
  for(unsigned int i = 0 ; i < 1 ; i++) {
    // unsigned int nz = 2*(i+1);
    unsigned int nz = 8;
    
    //double t_start = omp_get_wtime();
    mm_propagation.computeLawsonRK(nz);
  
    //double t_stop = omp_get_wtime();
    //double cpu_time = (t_stop - t_start);
  
    MultipleComplexArray phi_out(2,ComplexArray(nt));
    phi_out = mm_propagation.getResult();
    
    char filename[100];
    sprintf(filename,"results_GRIN2/cpp_nz%d.mat",nz);
    writeMatFile(phi_out,0.,filename);
    // writeMatFile(phi_out,cpu_time,"results_GRIN2/cpp_nz32768.mat");
  }
  
  return 0;
}

// ---------------------------------------------------------------------------- //

ComplexArray build_soliton(double t_final, unsigned int nt, double t0, double c, double q)
{
  ComplexArray phi(nt);
  
  complex<double> i(0.,1.);

  double h = (double) t_final/nt;
  double t = -t_final/2;

  for(unsigned int j = 0 ; j < nt ; j++) {
    phi[j] = 0.5 * sqrt(q/2) * (1./cosh(q/4 * (t-t0))) *  exp(i*c*((t-t0)/2));
    t += h;
  }
  return phi;
}

// ---------------------------------------------------------------------------- //
