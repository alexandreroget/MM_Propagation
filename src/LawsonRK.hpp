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

  void initializeLawson(const MultipleComplexArrays phi_in, const double delta_z);
  MultipleComplexArrays compute(const double delta_z, const unsigned int nz, const double z_stop);
   
private:
  void computeDifferentialOperator(const vector<doubleArray> beta);
  int factorial(int a);

private:
  unsigned int M; // Number of modes
  MultipleComplexArrays psi; // Lawson operator
 
  unsigned int nt;
  double t_final;

  unsigned int s;
  vector<double> c;

  MultipleComplexArrays L;
  
  vector<MultipleComplexArrays> E;
  MultipleComplexArrays E1;

  RungeKutta RK;
  Nonlinearity* N;
};

#endif
