#ifndef MULTIMODE_PROPAGATION_HPP
#define MULTIMODE_PROPAGATION_HPP

#include "LawsonRK.hpp"

typedef vector<double> doubleArray;
typedef vector<ComplexArray> MultipleComplexArrays;

struct SimulationParameters
{
  unsigned int n_modes;
  double fiber_length;
  unsigned int max_dispersion_order;
  vector<doubleArray> dispersion_coefficients;
  vector<Sparse3DMatrix> coupling_coefficients;
  
  unsigned int nt;
  double time_window;
  double pulse_width;
  MultipleComplexArrays signal;
  
  unsigned int n_steps;
  unsigned int order_RK;
  double nonlinearity_const;
  double raman_proportion;
  double raman_parameters[2];
};


class MultimodePropagation
{
public:
  MultimodePropagation(struct SimulationParameters initial_field);
  ~MultimodePropagation() { delete solver; }
  
  void computeLawsonRK(unsigned int n);  
  void set_nz(unsigned int n_steps);
  
  MultipleComplexArrays getResult() const;
  void saveResults(const MultipleComplexArrays phi, const string filename) const;
  
private:
  double get_max(const vector<double> beta) const;
  double get_max(const vector<Sparse3DMatrix> Q) const;

  void set_beta_tilde(vector<doubleArray>& beta_tilde, const vector<doubleArray> beta, const double T0, const double max_beta2);
  void set_Q_tilde(vector<Sparse3DMatrix>& Q_tilde, const vector<Sparse3DMatrix> Q, const double max_Q);
  
private:
  unsigned int M;

  unsigned int nt;
  double time_window;
  
  unsigned int nz;
  double z_final;

  double conversion_factor;
  
  MultipleComplexArrays phi_in;
  MultipleComplexArrays phi_out;
  
  LawsonRK* solver;
};

#endif
