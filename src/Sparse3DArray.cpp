#include "Sparse3DArray.hpp"

Sparse3DArray::Sparse3DArray() : num_nzv(0) {}


NonZeroValue Sparse3DArray::getNonZeroValue(const unsigned int i) const {
  return Q[i];
}
    

unsigned int Sparse3DArray::getNumberOfNonZeroValues() const {
  return num_nzv;
}


void Sparse3DArray::addNonZeroValue(const NonZeroValue nzv) {  
  Q.push_back(nzv);
  num_nzv++;
}
  

void Sparse3DArray::resetValue(const unsigned int i, const double Q_klm) {
  Q[i].value = Q_klm;
}
  

Sparse3DArray& Sparse3DArray::operator=(const Sparse3DArray& a) {
  Q.clear();
  num_nzv = a.num_nzv;
  
  for(unsigned int i = 0 ; i < num_nzv ; i++) {
    Q.push_back(a.Q[i]);
  }
  return *this;
}
