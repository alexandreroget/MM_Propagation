#ifndef COMPLEX_ARRAYS_CONTAINER_GPU_CUH_INCLUDED
#define COMPLEX_ARRAYS_CONTAINER_GPU_CUH_INCLUDED

// class ComplexArraysContainer;

#include <cuda_runtime.h>
#include <cuComplex.h>

#include "ComplexArraysContainer.hpp"

typedef ComplexArraysContainer ComplexArraysContainerCPU;


__global__ void copyManyComplexArraysKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);
__global__ void copySingleComplexArrayKernel(cuDoubleComplex* A, const unsigned int rowIndexA, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int rowIndexB, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);

__global__ void addManyComplexArraysKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);
__global__ void addSingleComplexArrayKernel(cuDoubleComplex* A, const unsigned int rowIndexA, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int rowIndexB, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);

__global__ void multiplyManyComplexArraysKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);
__global__ void multiplySingleComplexArrayKernel(cuDoubleComplex* A, const unsigned int rowIndexA, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int rowIndexB, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);
__global__ void multiplyManyComplexArraysByScalarKernel(cuDoubleComplex* A, const unsigned int columnIndex, const cuDoubleComplex b, const unsigned int arraySize, const unsigned int numArrays);

__global__ void computeComplexConjugateKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays);


class ComplexArraysContainerGPU
{
private:
  unsigned int nRows;
  unsigned int nCols;
  
  unsigned int n;
  
  cuDoubleComplex* data;
  
public:
  ComplexArraysContainerGPU() {}
  ComplexArraysContainerGPU(unsigned int arraySize, unsigned int numRows, unsigned int numColumns = 1);
  ComplexArraysContainerGPU(unsigned int arraySize, cuDoubleComplex value, unsigned int numRows, unsigned int numColumns = 1);
  ~ComplexArraysContainerGPU();
  
  void getDataFromCPU(ComplexArraysContainerCPU& source, const unsigned int columnIndex = 0);
  void sendDataToCPU(ComplexArraysContainerCPU& target, const unsigned int columnIndex = 0) const;

  cuDoubleComplex* operator () () const;  
  operator cuDoubleComplex* () const;
};


#endif
