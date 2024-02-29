#ifndef MULTIMODE_PROPAGATION_GPU_CUH
#define MULTIMODE_PROPAGATION_GPU_CUH

#include "MultimodePropagation.hpp"
#include "ComplexArraysContainerGPU.cuh"

#include <cufft.h>


class MultimodePropagationGPU : public MultimodePropagation {
public:
  MultimodePropagationGPU(const struct SimulationParameters& in);
  ~MultimodePropagationGPU();
  
  void computeLawsonRK(const unsigned int nz);
  
  ComplexArraysContainerCPU getResult();
  
protected:
  void setGPUThreadConfiguration();

  void computeExponentialOperatorsMatrix(ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& E1_GPU);
  
  void applyRungeKuttaMethod(ComplexArraysContainerGPU& psi_GPU, const ComplexArraysContainerGPU& E_GPU, double h);
  virtual void computeNonlinearity(ComplexArraysContainerGPU& N, const ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& psi_i, const unsigned int i) = 0;
  
  virtual void setFFTPlans() = 0;
  void setFFTPlansMany();
    
  void computeManyFFT(ComplexArraysContainerGPU& u);
  void computeManyIFFT(ComplexArraysContainerGPU& u);
    
  void copyManyComplexArrays(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);
  void copySingleComplexArray(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int rowIndexA, const unsigned int rowIndexB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);
  
  void addManyComplexArrays(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);
  void addSingleComplexArray(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int rowIndexA, const unsigned int rowIndexB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);
  
  void multiplyManyComplexArrays(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);
  void multiplySingleComplexArray(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int rowIndexA, const unsigned int rowIndexB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);
  void multiplyManyComplexArraysByScalar(ComplexArraysContainerGPU& arraysContainerA, const cuDoubleComplex complexScalar, const unsigned int columnIndexA = 0);
  
  void computeComplexConjugate(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA = 0, const unsigned int columnIndexB = 0);

protected:
  int n_threads;
  int n_blocks;
  
  cufftHandle forward_plan_many;
  cufftHandle inverse_plan_many;
};


class MultimodePropagationGPU_RamanOFF : public MultimodePropagationGPU {
public:
  MultimodePropagationGPU_RamanOFF(const struct SimulationParameters& in);
  ~MultimodePropagationGPU_RamanOFF() {}
  
private:
  void setFFTPlans();

  void computeNonlinearity(ComplexArraysContainerGPU& N, const ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& psi_i, const unsigned int i);
};


class MultimodePropagationGPU_RamanON : public MultimodePropagationGPU {
public:
  MultimodePropagationGPU_RamanON(const struct SimulationParameters& in);
  ~MultimodePropagationGPU_RamanON() {}
  
private:
  void setFFTPlans();

  void computeNonlinearity(ComplexArraysContainerGPU& N, const ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& psi_i, const unsigned int i);
  void convolveWithRamanResponse(ComplexArraysContainerGPU& V_Raman, const unsigned int index);
  
  void computeSingleFFT(ComplexArraysContainerGPU& u);
  void computeSingleIFFT(ComplexArraysContainerGPU& u);
  
private:
  double fR;
  ComplexArraysContainerGPU hR_fft;
  
  cufftHandle plan_1d;
};

#endif
