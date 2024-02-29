#ifndef MULTIMODE_PROPAGATION_HPP
#define MULTIMODE_PROPAGATION_HPP

#include <iostream>
#include <fstream>
#include <algorithm>

#include "ComplexArraysContainer.hpp"

#include "RungeKutta.hpp"
#include "Sparse3DArray.hpp"

typedef std::vector<double> doubleArray;


struct SimulationParameters
{
  unsigned int n_modes;
  double fiber_length;
  std::vector<doubleArray> dispersion_coefficients;
  std::vector<Sparse3DArray> coupling_coefficients;
  
  unsigned int nt;
  double time_window;
  double pulse_width;
  ComplexArraysContainer initial_fields;
  
  unsigned int n_steps;
  unsigned int method_order;
  double nonlinearity_const;
  double raman_proportion;
  doubleArray raman_response;
};


class MultimodePropagation { 
public:
  MultimodePropagation(const struct SimulationParameters& in);
  ~MultimodePropagation() {}
  
  virtual void computeLawsonRK(const unsigned int nz) = 0;
  void computeLawsonRK();
  virtual ComplexArraysContainer getResult() = 0;
  
protected:
  void setRungeKuttaCoefficients(const unsigned int method_order);
  void setDimensionlessVariables(const struct SimulationParameters& in);
  void setDimensionlessInitialFields(ComplexArraysContainer initial_fields);
  
  double getMax(const std::vector<doubleArray>& beta) const;
  double getMax(const std::vector<Sparse3DArray>& Q) const;
  
  std::vector<doubleArray> setDimensionlessDispersionCoefficients(const std::vector<doubleArray>& beta, const double LD, const double max_beta2);
  std::vector<Sparse3DArray> setDimensionlessCouplingCoefficients(const std::vector<Sparse3DArray>& Q, const double max_Q);
  
  void initLawsonRK();
  
  void computeLinearOperators();
  void initExponentialOperatorsMatrix();

  void setUniquePairs();
  virtual void setFFTPlans() = 0;
  
  int factorial(const int a);

protected:
  unsigned int M;
  unsigned int nt;
  
  double t_final;
  double delta_t;
  
  double z;
  double z_final;
  double delta_z;

  std::vector<doubleArray> beta_tilde;
  std::vector<Sparse3DArray> Q_tilde; 
  
  double T0;
  double sign_gamma;
  
  double conversion_factor;

  struct RungeKuttaCoefficients RK;
  
  ComplexArraysContainer L;
  
  ComplexArraysContainer psi;
  ComplexArraysContainer E1;
  std::vector<ComplexArraysContainer> E;
  
  std::vector<unsigned int> unique_pairs;
};

#endif
