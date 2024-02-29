#ifndef RUNGE_KUTTA_HPP
#define RUNGE_KUTTA_HPP

#include <math.h>
#include <vector>


typedef std::vector<double> doubleArray;

struct RungeKuttaCoefficients {
  unsigned int s;

  std::vector<doubleArray> a;
  doubleArray b;
  doubleArray c;
};


void Butcher4(struct RungeKuttaCoefficients& RK);
void Butcher7(struct RungeKuttaCoefficients& RK, const double nu);
void Butcher7(struct RungeKuttaCoefficients& RK, const double lambda, const double mu);
void Butcher11(struct RungeKuttaCoefficients& RK);

#endif
