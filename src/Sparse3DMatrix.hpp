#ifndef SPARSE_3D_MATRIX_HPP_INCLUDED
#define SPARSE_3D_MATRIX_HPP_INCLUDED

#include <vector>

using namespace std;


struct NonZeroValue
{
  unsigned int index[3];
  double value;
};

class Sparse3DMatrix
{
public:
  Sparse3DMatrix();

  void addNonZeroValue(const unsigned int k, const unsigned int l, const unsigned int m, const double Q_klm);
  
  unsigned int getSize() const;
  void getNonZeroValue(const unsigned int i, unsigned int &k, unsigned int &l, unsigned int &m, double &Q_klm) const;
  double getNonZeroValue(const unsigned int i) const;
 
  Sparse3DMatrix& operator=(const Sparse3DMatrix&);

private:
  unsigned int n;
  vector<struct NonZeroValue> Q;
};

#endif
