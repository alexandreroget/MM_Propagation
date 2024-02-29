#ifndef MULTIMODE_PROPAGATION_CPU_HPP
#define MULTIMODE_PROPAGATION_CPU_HPP

#include "MultimodePropagation.hpp"
#include "DiscreteFourierTransform.hpp"


class MultimodePropagationCPU : public MultimodePropagation {
public:
  MultimodePropagationCPU(const struct SimulationParameters& in);
  ~MultimodePropagationCPU() {}
  
  void computeLawsonRK(const unsigned int nz);
  
  ComplexArraysContainer getResult();

protected:
  void computeExponentialOperatorsMatrix();

  void applyRungeKuttaMethod(const double h);
  virtual void computeNonlinearity(ComplexArraysContainer& N_i, ComplexArraysContainer& E_i, ComplexArraysContainer& psi_i) = 0;
  
  void setFFTPlans();

protected:
  Forward_DFT** fft;
  Inverse_DFT** ifft;
};


class MultimodePropagationCPU_RamanOFF : public MultimodePropagationCPU {
public:
  MultimodePropagationCPU_RamanOFF(const struct SimulationParameters& in);
  ~MultimodePropagationCPU_RamanOFF() {}
  
private:
  void computeNonlinearity(ComplexArraysContainer& N_i, ComplexArraysContainer& E_i, ComplexArraysContainer& psi_i);
};


class MultimodePropagationCPU_RamanON : public MultimodePropagationCPU {
public:
  MultimodePropagationCPU_RamanON(const struct SimulationParameters& in);
  ~MultimodePropagationCPU_RamanON() {}
    
private:
  void computeNonlinearity(ComplexArraysContainer& N_i, ComplexArraysContainer& E_i, ComplexArraysContainer& psi_i);
  ComplexArray convolveWithRamanResponse(const ComplexArray& V);

private:
  double fR;
  ComplexArray hR_fft;
};

#endif
