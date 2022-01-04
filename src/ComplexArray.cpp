#include "ComplexArray.hpp"

using namespace std;

typedef complex<double> Complex;


// Constructor
ComplexArray::ComplexArray(unsigned int size) : n(size)
{
  v = new Complex[n];
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] = Complex(0.,0.);
  }
}


ComplexArray::ComplexArray(unsigned int size, Complex a) : n(size)
{
  v = new Complex[n];
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] = a;
  }
}


// Copy constructor
ComplexArray::ComplexArray(ComplexArray const& a) : n(a.n)
{
  v = new Complex[n];
  for(unsigned int i = 0 ; i < n ; i++) 
  {
    v[i] = a[i];
  }
}


int ComplexArray::getSize() const
{
  return n;
}


ComplexArray& ComplexArray::operator = (ComplexArray const a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] = a.v[i];
  }
  return *this;
}

ComplexArray& ComplexArray::operator = (ComplexArray const* a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] = a->v[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator += (const ComplexArray& a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] += a.v[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator -= (const ComplexArray& a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] -= a.v[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator *= (const ComplexArray& a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] *= a.v[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator *= (const Complex a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] *= a;
  }
  return *this;
}


ComplexArray& ComplexArray::operator *= (double a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] *= a;
  }
  return *this;
}


ComplexArray& ComplexArray::operator /= (const ComplexArray& a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] /= a.v[i];
  }
  return *this;
}


ComplexArray& ComplexArray::operator /= (const Complex a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] /= a;
  }
  return *this;
}


ComplexArray& ComplexArray::operator /= (double a)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    v[i] /= a;
  }
  return *this;
}


ComplexArray operator + (ComplexArray const& a, ComplexArray const& b)
{
  ComplexArray copy(a);
  copy += b;
  return copy;
}


ComplexArray operator - (ComplexArray const& a, ComplexArray const& b)
{
  ComplexArray copy(a);
  copy -= b;
  return copy;
}


ComplexArray operator * (ComplexArray const& a, ComplexArray const& b)
{
  ComplexArray copy(a);
  copy *= b;
  return copy;
}


ComplexArray operator * (ComplexArray const& a, Complex const& b)
{
  ComplexArray copy(a);
  copy *= b;
  return copy;
}


ComplexArray operator * (ComplexArray const& a, double const& b)
{
  ComplexArray copy(a);
  copy *= b;
  return copy;
}


ComplexArray operator / (ComplexArray const& a, ComplexArray const& b)
{
  ComplexArray copy(a);
  copy /= b;
  return copy;
}


ComplexArray operator / (ComplexArray const& a, Complex const& b)
{
  ComplexArray copy(a);
  copy /= b;
  return copy;
}


ComplexArray operator / (ComplexArray const& a, double const& b)
{
  ComplexArray copy(a);
  copy /= b;
  return copy;
}


double sum(const vector<double> A)
{
  int n = A.size();

  double result = 0.;
  for(unsigned int i = 0 ; i < n ; i++) {
    result += A[i];  
  }
  return result;
}

double norm(const ComplexArray A)
{
  int n = A.getSize();

  double result = 0;
  for(unsigned int i = 0 ; i < n ; i++) {
    result += std::norm(A[i]);
  }
  return sqrt(result);
}


vector<double> abs(const ComplexArray A)
{
  int n = A.getSize();

  vector<double> result;
  for(unsigned int i = 0 ; i < n ; i++) {
    result.push_back(std::abs(A[i]));   
  }
  return result;
}


ComplexArray conj(const ComplexArray A)
{
  int n = A.getSize();

  ComplexArray A_conj(n);
  for(unsigned int i = 0 ; i < n ; i++) {
    A_conj[i] = std::conj(A[i]);
  }
  return A_conj;
}


ComplexArray exp(const ComplexArray A)
{
  int n = A.getSize();

  ComplexArray exp_A(n);
  for(unsigned int i = 0 ; i < n ; i++) {
    exp_A[i] = std::exp(A[i]);
  }
  return exp_A;
}
