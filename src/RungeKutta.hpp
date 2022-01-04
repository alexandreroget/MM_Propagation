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
  void Butcher6(const double upsilon);

  unsigned int get_order() const;
  void get_c(vector<double>& coeff) const;

  MultipleComplexArrays apply_method(const double h, const vector<MultipleComplexArrays> E, 
                                     const MultipleComplexArrays psi, Nonlinearity* f);
  
private:
  unsigned int M;
  unsigned int nt;

  unsigned int s;

  vector<vector<double>> a;
  vector<double> b;
  vector<double> c;
};

#endif
