#ifndef COMPLEX_ARRAYS_CONTAINER_HPP_INCLUDED
#define COMPLEX_ARRAYS_CONTAINER_HPP_INCLUDED

#include "ComplexArray.hpp"


class ComplexArraysContainer
{
private:
  unsigned int M;
  unsigned int n;
  
  std::vector<ComplexArray> arrays;

public:
  ComplexArraysContainer() {}
  ComplexArraysContainer(unsigned int arraySize, unsigned int numArray);
  ComplexArraysContainer(unsigned int arraySize, std::complex<double> value, unsigned int numArrays);
  
  unsigned int getNumberOfArrays() const;
  unsigned int getArraySize() const;
  std::complex<double> getValue(const unsigned int arrayIndex, const unsigned int i) const;
  
  ComplexArray& operator [] (unsigned int i);
  ComplexArray& operator () (unsigned int i);
  
  std::vector<ComplexArray> operator () () const;  
  operator std::vector<ComplexArray> () const;
};

#endif
