#include "RungeKutta.hpp"

typedef vector<ComplexArray> MultipleComplexArrays;


// Constructor
RungeKutta::RungeKutta(unsigned int M, unsigned int nt, unsigned int order) : M(M), nt(nt)
{
  switch(order) {
    case 6:
      s = 7;
      Butcher7(1.);
      break;
    default:
      s = 4;
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


void RungeKutta::Butcher7(const double nu)
{
  a.push_back({0.});
  
  a.push_back({nu});
  
  a.push_back({(4.*nu-1.)/(8.*nu), 
                1./(8.*nu)});
  
  a.push_back({(10.*nu-2.)/(27.*nu), 
                2./(27.*nu), 
                (8.*nu)/(27.*nu)});
  
  a.push_back({(56.-77.*nu+(8.-17.*nu)*sqrt(21.))/(392.*nu), 
               (-8.*(7.+sqrt(21.)))/(392.*nu), 
               (48.*(7.+sqrt(21.))*nu)/(392.*nu), 
               (-3.*(21.+sqrt(21.))*nu)/(392.*nu)});
  
  a.push_back({(-5.*(287.*nu-56.-(59.*nu-8.)*sqrt(21.)))/(1960.*nu), 
               (-40.*(7.-sqrt(21.)))/(1960.*nu), 
               (320.*sqrt(21.)*nu)/(1960.*nu), 
               (3.*(21.-121.*sqrt(21.))*nu)/(1960.*nu), 
               (392.*(6.-sqrt(21.))*nu)/(1960.*nu)});
  
  a.push_back({(15.*((30.*nu-8.)-7.*nu*sqrt(21.)))/(180.*nu), 
               120./(180.*nu), 
               (-40.*(5.+7.*sqrt(21.))*nu)/(180.*nu), 
               (63.*(2.+3.*sqrt(21.))*nu)/(180.*nu), 
               (-14.*(49.-9.*sqrt(21.))*nu)/(180.*nu), 
               (70.*(7.+sqrt(21.))*nu)/(180.*nu)});
  
  b.push_back(9/180.);
  b.push_back(0.);
  b.push_back(64/180.);
  b.push_back(0.);
  b.push_back(49/180.);
  b.push_back(49/180.);
  b.push_back(9/180.);

  c.push_back(0.);
  c.push_back(nu);
  c.push_back(1/2.);
  c.push_back(2/3.);
  c.push_back((7.+sqrt(21.))/14.);
  c.push_back((7.-sqrt(21.))/14.);
  c.push_back(1.);
}


void RungeKutta::Butcher7(const double lambda, const double mu)
{
  a.push_back({0.});
  
  a.push_back({mu});
  
  a.push_back({2/3.-2./(9.*mu),
               2./(9.*mu)});
  
  a.push_back({5/12.-1./(9*mu), 
               1./(9*mu), 
               -1/12.});
  
  a.push_back({17/16.-3./(8.*mu), 
               3./(8.*mu), 
               -3/16., 
               -3/8.});
  
  a.push_back({17/16.-3./(8*mu)+1./(4.*lambda), 
               3./(8.*mu), 
               -3/16.-3./(4.*lambda), 
               -3/8.-3./(2.*lambda), 
               2./lambda});
  
  a.push_back({-27/44.+3./(11.*mu), 
               -3./(11.*mu), 63/44., 
               18/11., 
               4.*((lambda-4.)/11.), 
               -(4.*lambda)/11.});
  
  b.push_back(1/120.);
  b.push_back(0.);
  b.push_back(27/40.);
  b.push_back(27/40.);
  b.push_back((lambda-8.)/15.);
  b.push_back(-lambda/15.);
  b.push_back(11/120.);
  
  c.push_back(0.);
  c.push_back(mu);
  c.push_back(2/3.);
  c.push_back(1/3.);
  c.push_back(1/2.);
  c.push_back(1/2.);
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


// Apply Runge Kutta method
double RungeKutta::apply_method(const double h, const double lambda, const double y)
{  
  vector<double> f;
  double result = y;
  
  for(unsigned int i = 0 ; i < s ; i++) {
    double y_ni = y;

    for(unsigned int j = 0 ; j < i ; j++) {
      if(a[i][j] != 0.) {
        y_ni += (f[j]*(a[i][j]*h));
      }
    }

    f.push_back(lambda*y_ni);

    if(b[i] != 0.) {
      result += f[i]*(b[i]*h);
    }
  }
  
  return result;
}


// Apply Runge Kutta method
double RungeKutta::apply_method(const double h, const double z_n, const double lambda, const double psi)
{  
  vector<double> f(s);
  
  double result = psi;
  
  for(unsigned int i = 0 ; i < s ; i++) {
    double psi_ni = psi;
    
    for(unsigned int j = 0 ; j < i ; j++) {
      if(a[i][j] != 0.) {
        psi_ni += (f[j]*(a[i][j]*h));
      }
    }

    f[i] = (lambda*cos(z_n+c[i]*h)*exp(lambda*sin(z_n+c[i]*h)))*psi_ni;

    // psi[p]_{n+1} += b_i*h*N_op_{n,i}
    if(b[i] != 0.) {
      result += f[i]*(b[i]*h);
    }
  }
  
  return result;
}
