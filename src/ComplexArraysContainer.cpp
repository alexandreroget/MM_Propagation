#include "ComplexArraysContainer.hpp"


ComplexArraysContainer::ComplexArraysContainer(unsigned int arraySize, unsigned int numArrays) : M(numArrays), n(arraySize) {
  arrays = std::vector<ComplexArray>(M, ComplexArray(n));
}


ComplexArraysContainer::ComplexArraysContainer(unsigned int arraySize, std::complex<double> value, unsigned int numArrays) : M(numArrays), n(arraySize) {
  arrays = std::vector<ComplexArray>(M, ComplexArray(n, value));
}


unsigned int ComplexArraysContainer::getNumberOfArrays() const {
  return M;
}


unsigned int ComplexArraysContainer::getArraySize() const {
  return n;
}


std::complex<double> ComplexArraysContainer::getValue(const unsigned int arrayIndex, const unsigned int i) const {
  return arrays[arrayIndex][i];
}


ComplexArray& ComplexArraysContainer::operator [] (unsigned int i) { return arrays[i]; }
ComplexArray& ComplexArraysContainer::operator () (unsigned int i) { return arrays[i]; }

std::vector<ComplexArray> ComplexArraysContainer::operator () () const { return arrays; }
ComplexArraysContainer::operator std::vector<ComplexArray> () const { return arrays; }
