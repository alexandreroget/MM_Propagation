#ifndef RUNGEKUTTA_HPP_INCLUDED
#define RUNGEKUTTA_HPP_INCLUDED

#include "Nonlinearity.hpp"

using namespace std;

typedef vector<ComplexArray> MultipleComplexArrays;


class RungeKutta
{
public:
  RungeKutta(unsigned int M, unsigned int nt, unsigned int order);
  
  // void Butcher2(const double alpha);
  void Butcher4();
  void Butcher7(const double nu);
  void Butcher7(const double lambda, const double mu);

  unsigned int get_order() const;
  void get_c(vector<double>& coeff) const;

  MultipleComplexArrays apply_method(const double h, const vector<MultipleComplexArrays> E, 
                                     const MultipleComplexArrays psi, Nonlinearity* f);
                                     
  double apply_method(const double h, const double lambda, const double y);
  double apply_method(const double h, const double z_n, const double lambda, const double psi);
  
private:
  unsigned int M;
  unsigned int nt;

  unsigned int s;

  vector<vector<double>> a;
  vector<double> b;
  vector<double> c;
};

#endif
