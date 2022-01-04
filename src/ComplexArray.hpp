#ifndef COMPLEX_ARRAY_HPP_INCLUDED
#define COMPLEX_ARRAY_HPP_INCLUDED

#include <iostream>
#include <math.h>

#include <vector>
#include <complex>

using namespace std;

typedef complex<double> Complex;


class ComplexArray
{
public:
  ComplexArray(unsigned int size);
  ComplexArray(unsigned int size, Complex a);
  ComplexArray(Complex a);
  ComplexArray(ComplexArray const& a);
  
  ~ComplexArray() { delete v; }

  int getSize() const;

  Complex& operator [] (unsigned int i) const { return v[i]; }
  Complex& operator () (unsigned int i) const { return v[i]; }  
  Complex* operator () () const { return v; }

  operator Complex* () const { return v; }

  ComplexArray& operator = (ComplexArray const);
  ComplexArray& operator = (ComplexArray const*);
  
  ComplexArray& operator += (const ComplexArray&);
  ComplexArray& operator -= (const ComplexArray&);
  ComplexArray& operator *= (const ComplexArray&);
  ComplexArray& operator *= (const Complex);
  ComplexArray& operator *= (const double);
  ComplexArray& operator /= (const ComplexArray&);
  ComplexArray& operator /= (const Complex);
  ComplexArray& operator /= (const double);

private:
  int n;
  Complex* v;
};

ComplexArray operator + (ComplexArray const&, ComplexArray const&);
ComplexArray operator - (ComplexArray const&, ComplexArray const&);
ComplexArray operator * (ComplexArray const&, ComplexArray const&);
ComplexArray operator * (ComplexArray const&, double const&);
ComplexArray operator * (ComplexArray const&, Complex const&);
ComplexArray operator / (ComplexArray const&, double const&);
ComplexArray operator / (ComplexArray const&, Complex const&);
ComplexArray operator / (ComplexArray const&, double const&);

double sum(const vector<double>);
double norm(const ComplexArray);
vector<double> abs(const ComplexArray);
ComplexArray conj(const ComplexArray);
ComplexArray exp(const ComplexArray);

#endif
