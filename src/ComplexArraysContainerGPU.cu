#include "ComplexArraysContainerGPU.cuh"


__global__ void copyManyComplexArraysKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y;

  int indexA = (columnIndexA * (arraySize * numArrays)) + (idy * arraySize) + idx;
  int indexB = (columnIndexB * (arraySize * numArrays)) + (idy * arraySize) + idx;
  
  A[indexA] = B[indexB];
}


__global__ void copySingleComplexArrayKernel(cuDoubleComplex* A, const unsigned int rowIndexA, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int rowIndexB, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  unsigned int indexA = (columnIndexA * (arraySize * numArrays)) + (rowIndexA * arraySize) + idx;
  unsigned int indexB = (columnIndexB * (arraySize * numArrays)) + (rowIndexB * arraySize) + idx;
  
  A[indexA] = B[indexB];
}


__global__ void addManyComplexArraysKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y;

  int indexA = (columnIndexA * (arraySize * numArrays)) + (idy * arraySize) + idx;
  int indexB = (columnIndexB * (arraySize * numArrays)) + (idy * arraySize) + idx;
  
  A[indexA] = cuCadd(A[indexA], B[indexB]);
}


__global__ void addSingleComplexArrayKernel(cuDoubleComplex* A, const unsigned int rowIndexA, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int rowIndexB, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  unsigned int indexA = (columnIndexA * (arraySize * numArrays)) + (rowIndexA * arraySize) + idx;
  unsigned int indexB = (columnIndexB * (arraySize * numArrays)) + (rowIndexB * arraySize) + idx;
  
  A[indexA] = cuCadd(A[indexA], B[indexB]);
}


__global__ void multiplyManyComplexArraysKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y;
  
  int indexA = (columnIndexA * (arraySize * numArrays)) + (idy * arraySize) + idx;
  int indexB = (columnIndexB * (arraySize * numArrays)) + (idy * arraySize) + idx;
  
  A[indexA] = cuCmul(A[indexA], B[indexB]);
}


__global__ void multiplySingleComplexArrayKernel(cuDoubleComplex* A, const unsigned int rowIndexA, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int rowIndexB, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  
  unsigned int indexA = (columnIndexA * (arraySize * numArrays)) + (rowIndexA * arraySize) + idx;
  unsigned int indexB = (columnIndexB * (arraySize * numArrays)) + (rowIndexB * arraySize) + idx;
  
  A[indexA] = cuCmul(A[indexA], B[indexB]);
}


__global__ void multiplyManyComplexArraysByScalarKernel(cuDoubleComplex* A, const unsigned int columnIndex, const cuDoubleComplex b, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y;

  int index = (columnIndex * (arraySize * numArrays)) + (idy * arraySize) + idx;
  
  A[index] = cuCmul(A[index], b);
}


__global__ void computeComplexConjugateKernel(cuDoubleComplex* A, const unsigned int columnIndexA, const cuDoubleComplex* B, const unsigned int columnIndexB, const unsigned int arraySize, const unsigned int numArrays) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int idy = blockIdx.y * blockDim.y;

  int indexA = (columnIndexA * (arraySize * numArrays)) + (idy * arraySize) + idx;
  int indexB = (columnIndexB * (arraySize * numArrays)) + (idy * arraySize) + idx;
  
  A[indexA] = cuConj(B[indexB]);
}


ComplexArraysContainerGPU::ComplexArraysContainerGPU(unsigned int arraySize, unsigned int numRows, unsigned int numColumns) :
n(arraySize), nRows(numRows), nCols(numColumns) {
  unsigned int size = n * nRows * nCols;
  cudaMalloc((void**)&data, sizeof(cuDoubleComplex) * size);
}


ComplexArraysContainerGPU::ComplexArraysContainerGPU(unsigned int arraySize, cuDoubleComplex value, unsigned int numRows, unsigned int numColumns) :
n(arraySize), nRows(numRows), nCols(numColumns) {
  unsigned int size = n * nRows * nCols;
  cudaMalloc((void**)&data, sizeof(cuDoubleComplex) * size);
    
  cuDoubleComplex* host = new cuDoubleComplex[size];
  for(unsigned int i = 0 ; i < size ; i++) {
    host[i] = value;
  }
  
  cudaMemcpy(data, host, sizeof(cuDoubleComplex) * size, cudaMemcpyHostToDevice);
  
  delete[] host;
}


ComplexArraysContainerGPU::~ComplexArraysContainerGPU() {
  cudaFree(data);
}


void ComplexArraysContainerGPU::getDataFromCPU(ComplexArraysContainerCPU& source, const unsigned int columnIndex) {
  unsigned int size = n * nRows * nCols;
  cuDoubleComplex* host = new cuDoubleComplex[size];
  
  cudaMemcpy(host, *this, sizeof(cuDoubleComplex) * size, cudaMemcpyDeviceToHost);
  
  for(unsigned int p = 0 ; p < nRows ; p++) {
    unsigned int hostIndex = (columnIndex * nRows * nCols) + (p * n);
    for(unsigned int i = 0 ; i < n ; i++) {
      host[hostIndex] = make_cuDoubleComplex(source[p][i].real(), source[p][i].imag());
      hostIndex++;
    }
  }
  
  cudaMemcpy(*this, host, sizeof(cuDoubleComplex) * size, cudaMemcpyHostToDevice);

  delete[] host;
}


void ComplexArraysContainerGPU::sendDataToCPU(ComplexArraysContainerCPU& target, const unsigned int columnIndex) const {
  unsigned int size = n * nRows * nCols;
  cuDoubleComplex* host = new cuDoubleComplex[size];
    
  cudaMemcpy(host, *this, sizeof(cuDoubleComplex) * size, cudaMemcpyDeviceToHost);
  
  for(unsigned int p = 0 ; p < nRows ; p++) {
    unsigned int hostIndex = (columnIndex * nRows * nCols) + (p * n);
    for(unsigned int i = 0 ; i < n ; i++) {
      target[p][i] = std::complex<double>(cuCreal(host[hostIndex]), cuCimag(host[hostIndex]));
      hostIndex++;
    }
  }
  
  delete[] host;
}


cuDoubleComplex* ComplexArraysContainerGPU::operator () () const { 
  return data; 
}


ComplexArraysContainerGPU::operator cuDoubleComplex* () const { 
  return data; 
}

