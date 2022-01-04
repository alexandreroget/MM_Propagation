#ifndef LAWSONRK_HPP_INCLUDED
#define LAWSONRK_HPP_INCLUDED

#include <fstream>

#include "RungeKutta.hpp"

using namespace std;

typedef vector<double> doubleArray;
typedef vector<ComplexArray> MultipleComplexArrays;


class LawsonRK
{
public:
  LawsonRK(unsigned int M, double pulse_width, unsigned int nt, double t_final,
  unsigned int order_RK, double nonlinearity_const, vector<doubleArray> dispersion_coeff,
  vector<Sparse3DMatrix> coupling_coeff, double raman_proportion, double raman_parameters[2]);
  
  ~LawsonRK();

  void compute(const MultipleComplexArrays phi_in, MultipleComplexArrays& phi_out, const double z_final, const double delta_z);
   
private:
  void compute_L(const vector<doubleArray> beta);

  void compute_E1(const double h);
  void initialize_E(const double h);
  void compute_E();

  int factorial(int a);

private:
  unsigned int M; // Number of modes
  MultipleComplexArrays psi; // Lawson operator
 
  unsigned int nt;
  double t_final;

  unsigned int order_RK;
  vector<double> c;

  MultipleComplexArrays L;

  vector<MultipleComplexArrays> E;
  MultipleComplexArrays E1;

  RungeKutta RK;
  Nonlinearity* N;
};

#endif
