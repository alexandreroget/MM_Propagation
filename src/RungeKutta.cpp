#include "RungeKutta.hpp"

typedef vector<ComplexArray> MultipleComplexArrays;


// Constructor
RungeKutta::RungeKutta(unsigned int M, unsigned int nt, unsigned int order) : M(M), nt(nt), s(order)
{
  switch(s) {
    case 6:
      Butcher6(1.);
      break;
    default:
      Butcher4();
      break;
  }
}

/*
void RungeKutta::Butcher2(const double alpha)
{
  a.push_back({0.});
  a.push_back({alpha});
  
  b.push_back(1.-1./(2.*alpha));
  b.push_back(1./(2.*alpha));

  c.push_back(0.);
  c.push_back(alpha);
}
*/

void RungeKutta::Butcher4()
{
  a.push_back({0.});
  a.push_back({1/2.});
  a.push_back({0., 1/2.});
  a.push_back({0., 0., 1.});

  b.push_back(1/6.);
  b.push_back(1/3.);
  b.push_back(1/3.);
  b.push_back(1/6.);
    
  c.push_back(0.);
  c.push_back(1/2.);
  c.push_back(1/2.);
  c.push_back(1.);
}


void RungeKutta::Butcher6(const double upsilon)
{
  a.push_back({0.});
  a.push_back({upsilon});
  a.push_back({(4.*upsilon-1.)/(8.*upsilon), 1./(8.*upsilon)});
  a.push_back({(10.*upsilon-2.)/(27.*upsilon), 2./(27.*upsilon), (8.*upsilon)/(27.*upsilon)});
  a.push_back({(56.-77.*upsilon+(8.-17.*upsilon)*sqrt(21.))/(392.*upsilon), (-8.*(7.+sqrt(21.))*upsilon)/(392.*upsilon), (48.*(7.+sqrt(21.))*upsilon)/(392.*upsilon), (-3.*(21.+sqrt(21.))*upsilon)/(392.*upsilon)});
  a.push_back({(-5.*(287.*upsilon-56.-(59.*upsilon-8.)*sqrt(21.)))/(1960.*upsilon), (-40.*(7.-sqrt(21.)))/(1960.*upsilon), (320.*sqrt(21.)*upsilon)/(1960.*upsilon), (3.*(21.-121.*sqrt(21.))*upsilon)/(1960.*upsilon), (392.*(6.-sqrt(21.))*upsilon)/(1960.*upsilon)});
  a.push_back({(15.*(30.*upsilon-8.)-7.*upsilon*sqrt(21.))/(180.*upsilon), 120./(180.*upsilon), (-40.*(5.+7.*sqrt(21.))*upsilon)/(180.*upsilon), (63.*(2.+3.*sqrt(21.))*upsilon)/(180.*upsilon), (-14.*(49.-9.*sqrt(21.))*upsilon)/(180.*upsilon), (70.*(7.+sqrt(2.1))*upsilon)/(180.*upsilon)});
  
  b.push_back(9./180.);
  b.push_back(0.);
  b.push_back(64./180.);
  b.push_back(0.);
  b.push_back(49./180.);
  b.push_back(49./180.);
  b.push_back(9./180.);
  
  c.push_back(0.);
  c.push_back(upsilon);
  c.push_back(1./2.);
  c.push_back(2./3.);
  c.push_back((7.+sqrt(21.))/14.);
  c.push_back((7.-sqrt(21.))/14.);
  c.push_back(1.);
}


unsigned int RungeKutta::get_order() const
{
  return s;
}


void RungeKutta::get_c(vector<double>& coeff) const
{
  for(unsigned int i = 0 ; i < s ; i++) {
    coeff.push_back(c[i]);
  }
}


// Apply Runge Kutta method
MultipleComplexArrays RungeKutta::apply_method(const double h, const vector<MultipleComplexArrays> E, 
                                               const MultipleComplexArrays psi, Nonlinearity* f)
{  
  vector<MultipleComplexArrays> N;
  
  // psi[p]_{n+1} = psi[p]
  MultipleComplexArrays result;
  for(unsigned int p = 0 ; p < M ; p++) {
    result.push_back(psi[p]);
  }
  
  for(unsigned int i = 0 ; i < s ; i++) {
    // psi_{n,i} = psi_n + h*sum_{j=0}^{i-1} [a_{i,j} * N_op_{n,j}]
    MultipleComplexArrays psi_ni;
    for(unsigned int p = 0 ; p < M ; p++) {
      // Initialize psi_{n,i} = psi_n
      psi_ni.push_back(psi[p]);
    
      // if (i == 0) then psi_{n,i} = psi_n
      // else psi_{n,i} = psi_n + h*sum_{j=0}^{i-1} [a_{i,j} * N_op_{n,j}]
      for(unsigned int j = 0 ; j < i ; j++) {
        if(a[i][j] != 0.) {
          psi_ni[p] += (N[j][p]*(a[i][j]*h));
        }
      }
    }

    // N(z_n+ci*h,psi_ni)
    N.push_back(f->compute(E[i],psi_ni));

    // psi[p]_{n+1} += b_i*h*N_op_{n,i}
    if(b[i] != 0.) {
      for(unsigned int p = 0 ; p < M ; p++) {
        result[p] += N[i][p]*(b[i]*h);
      }
    }
  }
  
  return result;
}
