#ifndef COMPLEX_ARRAY_HPP_INCLUDED
#define COMPLEX_ARRAY_HPP_INCLUDED

#include <math.h>
#include <vector>
#include <complex>


typedef std::complex<double> Complex;


class ComplexArray
{
public:
  ComplexArray();
  ComplexArray(unsigned int size);
  ComplexArray(unsigned int size, Complex a);
  ComplexArray(ComplexArray const& a);
  
  ~ComplexArray();

  unsigned int getSize() const;

  Complex& operator [] (unsigned int i) const;
  Complex& operator () (unsigned int i) const;
  Complex* operator () () const;
  
  operator Complex* () const;

  ComplexArray& operator = (ComplexArray const a);
  ComplexArray& operator = (ComplexArray const* a);
  ComplexArray& operator += (const ComplexArray& a);
  ComplexArray& operator += (const Complex a);
  ComplexArray& operator -= (const ComplexArray& a);
  ComplexArray& operator *= (const ComplexArray& a);
  ComplexArray& operator *= (const Complex a);
  ComplexArray& operator *= (const double a);
  ComplexArray& operator /= (const ComplexArray& a);
  ComplexArray& operator /= (const Complex a);
  ComplexArray& operator /= (const double a);

private:
  unsigned int n;
  Complex* data;
};

ComplexArray operator + (ComplexArray const& a, ComplexArray const& b);
ComplexArray operator - (ComplexArray const& a, ComplexArray const& b);
ComplexArray operator * (ComplexArray const& a, ComplexArray const& b);
ComplexArray operator * (ComplexArray const& a, Complex const& b);
ComplexArray operator * (ComplexArray const& a, double const& b);
ComplexArray operator / (ComplexArray const& a, ComplexArray const& b);
ComplexArray operator / (ComplexArray const& a, Complex const& b);
ComplexArray operator / (ComplexArray const& a, double const& b);

double sum(const std::vector<double> A);
double norm(const ComplexArray A);
std::vector<double> abs(const ComplexArray A);
ComplexArray conj(const ComplexArray A);
ComplexArray exp(const ComplexArray A);

#endif
