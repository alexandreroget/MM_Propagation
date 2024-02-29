#include "ComplexArray.hpp"

ComplexArray::ComplexArray() : n(0) {
  data = NULL;
}


ComplexArray::ComplexArray(unsigned int size) : n(size) {
  data = new Complex[n];
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] = Complex(0.,0.);
  }
}


ComplexArray::ComplexArray(unsigned int size, Complex a) : n(size) {
  data = new Complex[n];
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] = a;
  }
}


ComplexArray::ComplexArray(ComplexArray const& a) : n(a.n) {
  data = new Complex[n];
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] = a[i];
  }
}


ComplexArray::~ComplexArray() { delete data; }

unsigned int ComplexArray::getSize() const {
  return n;
}


Complex& ComplexArray::operator [] (unsigned int i) const { return data[i]; }
Complex& ComplexArray::operator () (unsigned int i) const { return data[i]; }  
Complex* ComplexArray::operator () () const { return data; }
ComplexArray::operator Complex* () const { return data; }


ComplexArray& ComplexArray::operator = (ComplexArray const a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] = a.data[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator = (ComplexArray const* a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] = a->data[i];
  }
  return *this;
}

  
ComplexArray& ComplexArray::operator += (const ComplexArray& a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] += a.data[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator += (const Complex a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] += a;
  }
  return *this;
}


ComplexArray& ComplexArray::operator -= (const ComplexArray& a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] -= a.data[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator *= (const ComplexArray& a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] *= a.data[i];
  }
  return *this;
}

  
ComplexArray& ComplexArray::operator *= (const Complex a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] *= a;
  }
  return *this;
}
  

ComplexArray& ComplexArray::operator *= (const double a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] *= a;
  }
  return *this;
}
  

ComplexArray& ComplexArray::operator /= (const ComplexArray& a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] /= a.data[i];
  }
  return *this;
}
  

ComplexArray& ComplexArray::operator /= (const Complex a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] /= a;
  }
  return *this;
}
  

ComplexArray& ComplexArray::operator /= (const double a) {
  for(unsigned int i = 0 ; i < n ; i++) {
    data[i] /= a;
  }
  return *this;
}


ComplexArray operator + (ComplexArray const& a, ComplexArray const& b) {
  ComplexArray copy(a);
  copy += b;
  return copy;
}

ComplexArray operator - (ComplexArray const& a, ComplexArray const& b) {
  ComplexArray copy(a);
  copy -= b;
  return copy;
}

ComplexArray operator * (ComplexArray const& a, ComplexArray const& b) {
  ComplexArray copy(a);
  copy *= b;
  return copy;
}

ComplexArray operator * (ComplexArray const& a, Complex const& b) {
  ComplexArray copy(a);
  copy *= b;
  return copy;
}

ComplexArray operator * (ComplexArray const& a, double const& b) {
  ComplexArray copy(a);
  copy *= b;
  return copy;
}

ComplexArray operator / (ComplexArray const& a, ComplexArray const& b) {
  ComplexArray copy(a);
  copy /= b;
  return copy;
}


ComplexArray operator / (ComplexArray const& a, Complex const& b) {
  ComplexArray copy(a);
  copy /= b;
  return copy;
}


ComplexArray operator / (ComplexArray const& a, double const& b) {
  ComplexArray copy(a);
  copy /= b;
  return copy;
}


double sum(const std::vector<double> A) {
  unsigned int n = A.size();

  double result = 0.;
  for(unsigned int i = 0 ; i < n ; i++) {
    result += A[i];  
  }
  return result;
}


double norm(const ComplexArray A) {
  unsigned int n = A.getSize();

  double result = 0;
  for(unsigned int i = 0 ; i < n ; i++) {
    result += std::norm(A[i]);
  }
  return sqrt(result);
}


std::vector<double> abs(const ComplexArray A) {
  unsigned int n = A.getSize();

  std::vector<double> result;
  for(unsigned int i = 0 ; i < n ; i++) {
    result.push_back(std::abs(A[i]));   
  }
  return result;
}


ComplexArray conj(const ComplexArray A) {
  unsigned int n = A.getSize();

  ComplexArray A_conj(n);
  for(unsigned int i = 0 ; i < n ; i++) {
    A_conj[i] = std::conj(A[i]);
  }
  return A_conj;
}


ComplexArray exp(const ComplexArray A) {
  unsigned int n = A.getSize();

  ComplexArray exp_A(n);
  for(unsigned int i = 0 ; i < n ; i++) {
    exp_A[i] = std::exp(A[i]);
  }
  return exp_A;
}

