#include "Sparse3DMatrix.hpp"


Sparse3DMatrix::Sparse3DMatrix() : n(0)
{
}


void Sparse3DMatrix::addNonZeroValue(const unsigned int k, const unsigned int l, const unsigned int m, const double Q_klm)
{
  struct NonZeroValue nzv;
  
  nzv.index[0] = k;
  nzv.index[1] = l;
  nzv.index[2] = m;
  nzv.value = Q_klm;
  
  Q.push_back(nzv);
  n++;
}


unsigned int Sparse3DMatrix::getSize() const
{
  return n;
}


void Sparse3DMatrix::getNonZeroValue(const unsigned int i, unsigned int &k, unsigned int &l, unsigned int &m, double &Q_klm) const
{
  k = Q[i].index[0];
  l = Q[i].index[1];
  m = Q[i].index[2];
  Q_klm = Q[i].value;
}


double Sparse3DMatrix::getNonZeroValue(const unsigned int i) const
{
  return Q[i].value;
}


Sparse3DMatrix& Sparse3DMatrix::operator=(const Sparse3DMatrix& a)
{
  Q.clear();
  n = a.n;
  
  for(unsigned int i = 0 ; i < n ; i++) {
    Q.push_back(a.Q[i]);
  }
  return *this;
}
