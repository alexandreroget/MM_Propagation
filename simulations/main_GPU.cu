#include <iostream>
#include <fstream>

#include "MultimodePropagation_GPU.cuh"

ComplexArray build_soliton(const double t_final, const unsigned int nt, const double t0, const double c, const double q);

// ---------------------------------------------------------------------------- //

int main() {
  unsigned int nt = (int) pow(2,15);
  double t_final = 160.; // in ps

  double t0 = -t_final/10;
  double c = 1.;

  // ----------- Set up input parameters ----------- //
  
  struct SimulationParameters in;
  
  in.n_modes = 2;
  
  in.dispersion_coefficients = std::vector<std::vector<double>>(2, std::vector<double>(3, 0.));
  in.dispersion_coefficients[0][2] = 18.938e-3;
  in.dispersion_coefficients[1][2] = 18.866e-3;
  
  double SR = 6.4341e12;
  
  Sparse3DArray Q0;
  struct NonZeroValue nzv;
  
  nzv = {0, 0, 0, SR};
  Q0.addNonZeroValue(nzv);
  nzv = {0, 0, 1, SR/3};
  Q0.addNonZeroValue(nzv);
  nzv = {0 , 1, 1, SR/2};
  Q0.addNonZeroValue(nzv);
  // nzv = {1, 1, 1, -SR/3};
  // Q0.addNonZeroValue(nzv);
  in.coupling_coefficients.push_back(Q0);
  
  Sparse3DArray Q1;
  // nzv = {0, 0, 0, -SR/3};
  // Q1.addNonZeroValue(nzv);
  nzv = {1, 0, 0, SR/5};
  Q1.addNonZeroValue(nzv);
  nzv = {1, 1, 0, SR/3};
  Q1.addNonZeroValue(nzv);
  nzv = {1, 1, 1, 0.8*SR};
  Q1.addNonZeroValue(nzv);
  in.coupling_coefficients.push_back(Q1);
  
  in.fiber_length = 16.;

  in.n_steps = 1024;
  
  in.method_order = 4;

  double lambda = 1030e-9; // central wavelength in m
  double f0 = 2.99792458e-4/lambda;
  double speed_of_light = 2.99792458e-4; // speed of light in m/ps
  double w0 = 2*M_PI*f0;
  double n2 = 2.3*1e-20; // m^2 W^-1
  in.nonlinearity_const = (n2*w0)/speed_of_light;
  
  in.raman_proportion = 0.18;
  
  in.raman_response = std::vector<double>(nt);
  double tau[2] = {12.2e-3, 32e-3};
  double delta = t_final/nt;
  double a = (1./(tau[0]*tau[0]) + 1./(tau[1]*tau[1])) * tau[0];
  for(unsigned int i = 0 ; i < nt ; i++) {
    double t = i * delta;
    double b = exp(-t/tau[1]) * sin(t/tau[0]);
    in.raman_response[i] = a*b;
  }

  double beta_2 = in.dispersion_coefficients[0][2];

  in.nt = nt;
  in.pulse_width = std::sqrt(beta_2)/(SR*in.nonlinearity_const);
  in.time_window = t_final;

  in.initial_fields = ComplexArraysContainer(nt, 2);
  in.initial_fields[0] = build_soliton(t_final/std::sqrt(beta_2/2), nt, t0/std::sqrt(beta_2/2), c*std::sqrt(2/beta_2), SR*in.nonlinearity_const);
  in.initial_fields[1] = build_soliton(t_final/std::sqrt(beta_2/2), nt, -t0/std::sqrt(beta_2/2), -c*std::sqrt(2/beta_2), 0.8*SR*in.nonlinearity_const);
  
  std::ofstream outputFile("GRIN2_GPU_RamanON.txt");
  
  if (!outputFile.is_open()) {
    std::cerr << "Erreur : Impossible d'ouvrir le fichier de sortie." << std::endl;
    return 1;
  }

  MultimodePropagation* mm_propagation = new MultimodePropagationGPU_RamanON(in);
  
  mm_propagation->computeLawsonRK(in.n_steps);
  
  ComplexArraysContainer Phi_out(nt, 2);
  Phi_out = mm_propagation->getResult();
  
  for(unsigned int i = 0 ; i < nt ; i++) {
    outputFile << std::abs(in.initial_fields[0][i]) << " ; " << std::abs(in.initial_fields[1][i]) << " ; " << std::abs(Phi_out[0][i]) << " ; " << std::abs(Phi_out[1][i]) << std::endl;
  }
  
  outputFile.close();
  
  delete mm_propagation;

  return 0;
}

// ---------------------------------------------------------------------------- //

ComplexArray build_soliton(const double t_final, const unsigned int nt, const double t0, const double c, const double q) {
  ComplexArray phi(nt);
  
  std::complex<double> i(0.,1.);

  double h = (double) t_final/nt;
  double t = -t_final/2;

  for(unsigned int j = 0 ; j < nt ; j++) {
    phi[j] = 0.5 * std::sqrt(q/2) * (1./std::cosh(q/4 * (t-t0))) *  std::exp(i*c*((t-t0)/2));
    t += h;
  }
  return phi;
}

