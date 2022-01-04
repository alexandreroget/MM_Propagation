#ifndef DISCRETE_FOURIER_TRANSFORM_INCLUDED
#define DISCRETE_FOURIER_TRANSFORM_INCLUDED

#include <fftw3.h>
#include <omp.h>

#include "ComplexArray.hpp"

typedef complex<double> Complex;


class DiscreteFourierTransform
{
public:
  DiscreteFourierTransform(unsigned int n_threads, unsigned int size, int direction);
  // virtual ~DiscreteFourierTransform() {}
  
  virtual void compute(Complex* in, Complex* out) = 0;
  
  void destroy_plan() { fftw_destroy_plan(plan); }
  
protected:
  unsigned int n;
  fftw_plan plan;
};


class Forward_DFT : public DiscreteFourierTransform
{
public:
  Forward_DFT(unsigned int n_threads, unsigned int size) : DiscreteFourierTransform(n_threads, size, FFTW_FORWARD) {}
  // ~Forward_DFT() { fftw_destroy_plan(plan); }
  void compute(Complex* in, Complex* out);
};


class Backward_DFT : public DiscreteFourierTransform
{
public:
  Backward_DFT(unsigned int n_threads, unsigned int size) : DiscreteFourierTransform(n_threads, size, FFTW_BACKWARD) {}
  // ~Backward_DFT() { fftw_destroy_plan(plan); }
  void compute(Complex* in, Complex* out);
};

#endif
