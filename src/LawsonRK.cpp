#include "LawsonRK.hpp"

using namespace std;

typedef vector<double> doubleArray;
typedef vector<ComplexArray> MultipleComplexArrays;


// Constructor
LawsonRK::LawsonRK(unsigned int M, double pulse_width, unsigned int nt, double t_final,
unsigned int order_RK, double nonlinearity_const, vector<doubleArray> dispersion_coeff,
vector<Sparse3DMatrix> coupling_coeff, double raman_proportion, double raman_parameters[2]) :
M(M), nt(nt), t_final(t_final), RK(M,nt,order_RK), psi(M,ComplexArray(nt)), L(M,ComplexArray(nt)), E1(M,ComplexArray(nt))
{
  switch(order_RK) {
    case 6:
      s = 7;
      break;
    default:
      s = 4;
      break;  
  }
  E = vector<MultipleComplexArrays>(s,MultipleComplexArrays(M,ComplexArray(nt)));

  computeDifferentialOperator(dispersion_coeff);

  RK.get_c(c);
  
  if(raman_proportion == 0.) {
    N = new Nonlinearity_withoutRaman(M,nt,t_final,nonlinearity_const,coupling_coeff,pulse_width);
  }
  else {
    N = new Nonlinearity_withRaman(M,nt,t_final,nonlinearity_const,coupling_coeff,pulse_width,raman_proportion,raman_parameters);
  }
}


LawsonRK::~LawsonRK()
{
  delete N;
}


void LawsonRK::initializeLawson(const MultipleComplexArrays phi_in, const double delta_z)
{
  for(unsigned int p = 0 ; p < M ; p++) {
    psi[p] = phi_in[p];

    E1[p] = exp(L[p]*delta_z);
    
    // compute E_{0,i}
    for(unsigned int i = 0 ; i < s ; i++) { 
      E[i][p] = exp(L[p]*(c[i]*delta_z));
    }
  }
}


MultipleComplexArrays LawsonRK::compute(const double delta_z, const unsigned int nz, const double z_stop)
{
  for(unsigned int i = 0 ; i < nz ; i++) {
    // compute psi_{n+1}
    psi = RK.apply_method(delta_z,E,psi,N);
    
    if(isnan(psi[0][0].real())) {
      cerr<<"Error: NaN field encountered."<<endl;
      exit(EXIT_FAILURE);
    }
    
    // E_{n+1,i} = E_{n,i} * E_{1,i}
    for(unsigned int i = 0 ; i < s ; i++) {
      for(unsigned int p = 0 ; p < M ; p++) {
        E[i][p] *= E1[p]; 
      }
    }
  }

  MultipleComplexArrays phi_out(M,ComplexArray(nt));

  unsigned int max_threads = omp_get_max_threads();
  Forward_DFT fft(max_threads,nt);
  Backward_DFT ifft(max_threads,nt);
  
  for(unsigned int p = 0 ; p < M ; p++) {
    ComplexArray psi_fft(nt);
    fft.compute(psi[p],psi_fft);
    
    psi_fft *= exp(L[p]*z_stop);
    ifft.compute(psi_fft,phi_out[p]);
  }

  fft.destroy_plan();
  ifft.destroy_plan();
  
  return phi_out;
}


void LawsonRK::computeDifferentialOperator(const vector<doubleArray> beta)
{
  complex<double> i(0.,1.);

  unsigned int n_beta = beta[0].size();
  
  vector<double> alpha(nt);
  for(unsigned int j = 0 ; j < nt ; j++) {
    int k;
    j<=nt/2 ? k=j : k=j-nt;
    alpha[j] = (2*M_PI/t_final) * k;
  }

  for(unsigned int p = 0 ; p < M ; p++) {
    for(unsigned int j = 0 ; j < nt ; j++) {
      // L = i*beta_0 - i*beta_1*alpha(k) - i*beta_2/2*alpha(k)^2
      L[p][j] = i * (beta[p][0] -beta[p][1]*alpha[j] -beta[p][2]/2*(alpha[j]*alpha[j]));
      if(n_beta >= 3) {
        for(unsigned int k = 3 ; k < n_beta ; k++) {
          L[p][j] += (i * (pow(-1,k) * beta[p][k]/factorial(k) * pow(alpha[j],k)));
        }
      }
    }
  }
}


int LawsonRK::factorial(int a)
{
  int result = 1;
  if(a > 0) {
    for(int i = 1 ; i <= a ; i++) {
      result *= i;
    }
  }
  return result;
}
