#ifndef SPARSE_3D_ARRAY_HPP_INCLUDED
#define SPARSE_3D_ARRAY_HPP_INCLUDED

#include <vector>


struct NonZeroValue
{
  unsigned int index[3];
  double value;
};


class Sparse3DArray
{
public:
  Sparse3DArray();

  NonZeroValue getNonZeroValue(const unsigned int i) const;
  unsigned int getNumberOfNonZeroValues() const;

  void addNonZeroValue(const NonZeroValue nzv);
  void resetValue(const unsigned int i, const double Q_klm);
  
  Sparse3DArray& operator=(const Sparse3DArray& a);

private:
  unsigned int num_nzv;
  std::vector<struct NonZeroValue> Q;
};

#endif
