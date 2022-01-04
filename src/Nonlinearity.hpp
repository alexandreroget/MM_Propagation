#ifndef NONLINEARITY_HPP_INCLUDED
#define NONLINEARITY_HPP_INCLUDED

#include <map>
#include <memory>

#include "Sparse3DMatrix.hpp"
#include "DiscreteFourierTransform.hpp"

using namespace std;

typedef vector<ComplexArray> MultipleComplexArrays;


class Nonlinearity
{
public:
  Nonlinearity(unsigned int M, unsigned int nt, double t_final,
  double nonlinearity_const, vector<Sparse3DMatrix> coupling_coeff, double pulse_width);
  
  virtual ~Nonlinearity() {}

  virtual MultipleComplexArrays compute(const MultipleComplexArrays E, const MultipleComplexArrays psi) = 0;

protected:
  unsigned int M;
  unsigned int nt;
  
  double t_final;
  double delta;
  
  // Nonlinarity index
  double gamma;
  int sign_gamma;

  // Pulse width
  double T0;
  
  // Coupling coefficients
  vector<Sparse3DMatrix> Q;

  vector<unsigned int> index_array;

  MultipleComplexArrays U;
  MultipleComplexArrays U_conj;
  MultipleComplexArrays N;

  Forward_DFT** fft;
  Backward_DFT** ifft;
};


class Nonlinearity_withoutRaman : public Nonlinearity
{
public:
  Nonlinearity_withoutRaman(unsigned int M, unsigned int nt, double t_final,
  double nonlinearity_const, vector<Sparse3DMatrix> coupling_coeff, double pulse_width);
  
  ~Nonlinearity_withoutRaman();

  MultipleComplexArrays compute(const MultipleComplexArrays E, const MultipleComplexArrays psi);
};


class Nonlinearity_withRaman : public Nonlinearity
{
public:
  Nonlinearity_withRaman(unsigned int n_modes, unsigned int n_time, double t_final,
  double nonlinearity_const, vector<Sparse3DMatrix> coupling_coeff,
  double pulse_width, double raman_proportion, double raman_parameters[2]);
  
  ~Nonlinearity_withRaman();

  MultipleComplexArrays compute(const MultipleComplexArrays E, const MultipleComplexArrays psi);

private:
  ComplexArray convolution_with_hR(ComplexArray& A);
  void compute_delayed_Raman_response(const double tau[2]);
  
  // Fractional contribution of the delayed Raman response  
  double fR;
  // FFT of the delayed Raman response
  ComplexArray hR_fft;
};

#endif
