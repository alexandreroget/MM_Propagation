#include "Nonlinearity.hpp"

typedef vector<ComplexArray> MultipleComplexArrays;


Nonlinearity::Nonlinearity(unsigned int M, unsigned int nt, double t_final, 
double nonlinearity_const, vector<Sparse3DMatrix> coupling_coeff, double pulse_width) :
M(M), nt(nt), t_final(t_final), gamma(nonlinearity_const), T0(pulse_width)
{
  delta = (double) t_final/nt;
  gamma >= 0 ? sign_gamma = 1 : sign_gamma = -1;

  unsigned int max_threads = omp_get_max_threads();
  unsigned int n_threads = max_threads/M;
  if(n_threads == 0) n_threads++;

  fft = new Forward_DFT*[M];
  ifft = new Backward_DFT*[M];
  
  for(unsigned int p = 0 ; p < M ; p++) {
    fft[p] = new Forward_DFT(n_threads,nt);
    ifft[p] = new Backward_DFT(n_threads,nt);
  
    U.push_back(ComplexArray(nt));
    U_conj.push_back(ComplexArray(nt));
    N.push_back(ComplexArray(nt));
  
    Q.push_back(Sparse3DMatrix());
    Q[p] = coupling_coeff[p];

    unsigned int n = Q[p].getSize();
    for(unsigned int i = 0 ; i < n ; i++) {
      unsigned int k, l, m;
      double Q_pklm;
      Q[p].getNonZeroValue(i,k,l,m,Q_pklm);
      unsigned int lm = M*l + m;
      
      bool new_lm = true;
      for(vector<unsigned int>::iterator it = index_array.begin() ; it != index_array.end() ; ++it) {
        if(*it == lm) {
          new_lm = false;
          break;
        }
      }
      if(new_lm) {
        index_array.push_back(lm);
      }
    }
  }
}


Nonlinearity_withoutRaman::~Nonlinearity_withoutRaman()
{
  for(unsigned int p = 0 ; p < M ; p++) {
    fft[p]->destroy_plan();
    delete[] fft[p];
    ifft[p]->destroy_plan();
    delete[] ifft[p];
  }
  delete[] fft;
  delete[] ifft;
}

Nonlinearity_withRaman::~Nonlinearity_withRaman()
{
  for(unsigned int p = 0 ; p < M ; p++) {
    fft[p]->destroy_plan();
    delete[] fft[p];
    ifft[p]->destroy_plan();
    delete[] ifft[p];
  }
  delete[] fft;
  delete[] ifft;
}

Nonlinearity_withoutRaman::Nonlinearity_withoutRaman(unsigned int M, unsigned int nt, double t_final,
double nonlinearity_const, vector<Sparse3DMatrix> coupling_coeff, double pulse_width) :
Nonlinearity(M,nt,t_final,nonlinearity_const,coupling_coeff,pulse_width)
{
}


MultipleComplexArrays Nonlinearity_withoutRaman::compute(MultipleComplexArrays E_i, MultipleComplexArrays psi)
{ 
  // Compute U = E_i * psi_ni
  #pragma omp parallel for shared(psi, E_i, U, U_conj, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {
    U[p] = psi[p];
    // Fast Fourier Transform of U
    ComplexArray U_fft(nt);
    fft[p]->compute(U[p],U_fft);
    // The product (psi_p_fft * E_i[p]) is made in the frequency domain
    U_fft *= E_i[p];
    // Inverse Fast Fourier Transform -> Back to the time domain
    ifft[p]->compute(U_fft,U[p]);
    // Conjugated form of T
    U_conj[p] = conj(U[p]);
  }
  
  unsigned int n = index_array.size();
  MultipleComplexArrays A;
  for(unsigned int j = 0 ; j < n ; j++) {
    unsigned int lm = index_array[j];
    unsigned int l = lm / M;
    unsigned int m = lm % M;

    A.push_back(U[l] * U_conj[m]);
  }
  
  // N = i * exp[(-z_n+ci*h)*Lp] * sum_{k,l,m} Q^{(p)}_{k,l,m} * (U[k] * U[l] * conj(U[m]))
  #pragma omp parallel for shared(Q, E_i, U, U_conj, N, index_array, A, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {  
    ComplexArray V(nt);
    ComplexArray sigma(nt);

    unsigned int n_couplings = Q[p].getSize();
    // sigma^{(p)} = sum_{k,l,m} Q^{(p)}_{k,l,m} * (U[k] * U[l] * conj(U[m]))
    for(unsigned int j = 0 ; j < n_couplings ; j++) {
      unsigned int k, l, m;
      double Q_pklm;
      Q[p].getNonZeroValue(j,k,l,m,Q_pklm);

      unsigned int lm = l*M + m;
      unsigned int i = 0;
      for(vector<unsigned int>::iterator it = index_array.begin() ; it != index_array.end() ; ++it) {
        if(lm == *it) {
          break;
        }
        i++;
      }

      V = A[i] * Q_pklm;
      V *= U[k];
      
      sigma += V;
    }
    
    // Fast Fourier Transform of the sum
    ComplexArray sigma_fft(nt);
    fft[p]->compute(sigma,sigma_fft);
    // The product (sigma * conj(E_i)) is made in the frequency domain
    // Then we go back to the time domain with an Inverse Fast Fourier Transform
    // Finally we multiply the result with i
    sigma_fft *= conj(E_i[p]);
    ifft[p]->compute(sigma_fft,N[p]);
    N[p] *= complex<double>(0.,sign_gamma);
  }

  return N;
}


Nonlinearity_withRaman::Nonlinearity_withRaman(unsigned int M, unsigned int nt, double t_final,
double nonlinearity_const, vector<Sparse3DMatrix> coupling_coeff, double pulse_width, 
double raman_proportion, double raman_parameters[2]) :
Nonlinearity(M,nt,t_final,nonlinearity_const,coupling_coeff,pulse_width), 
fR(raman_proportion), hR_fft(nt)
{
  compute_delayed_Raman_response(raman_parameters);
}

MultipleComplexArrays Nonlinearity_withRaman::compute(MultipleComplexArrays E_i, MultipleComplexArrays psi)
{ 
  // Compute U = E_i * psi_ni
  #pragma omp parallel for shared(psi, E_i, U, U_conj, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {
    U[p] = psi[p];
    // Fast Fourier Transform of U
    ComplexArray U_fft(nt);
    fft[p]->compute(U[p],U_fft);
    // The product (psi_p_fft * E_i[p]) is made in the frequency domain
    U_fft *= E_i[p];
    // Inverse Fast Fourier Transform -> Back to the time domain
    ifft[p]->compute(U_fft,U[p]);
    // Conjugated form of T
    U_conj[p] = conj(U[p]);
  }
  
  unsigned int n = index_array.size();
  MultipleComplexArrays A;
  MultipleComplexArrays hR_conv_A;
  for(unsigned int j = 0 ; j < n ; j++) {
    unsigned int lm = index_array[j];
    unsigned int l = lm / M;
    unsigned int m = lm % M;

    A.push_back(U[l] * U_conj[m]);
    hR_conv_A.push_back(convolution_with_hR(A[j]));
  }
  
  // N = i * exp[(-z_n+ci*h)*Lp] * sum_{k,l,m} Q^{(p)}_{k,l,m} * (Kerr^{(p)} + Raman^{(p)})
  #pragma omp parallel for shared(Q, E_i, U, U_conj, hR_fft, N, index_array, A, hR_conv_A, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {  
    ComplexArray V(nt);
    ComplexArray sigma(nt);

    unsigned int n_couplings = Q[p].getSize();
    // sigma^{(p)} = sum_{k,l,m} Q^{(p)}_{k,l,m} * (Kerr^{(p)} + Raman^{(p)})
    for(unsigned int j = 0 ; j < n_couplings ; j++) {
      unsigned int k, l, m;
      double Q_pklm;
      Q[p].getNonZeroValue(j,k,l,m,Q_pklm);

      unsigned int lm = l*M + m;
      unsigned int i = 0;
      for(vector<unsigned int>::iterator it = index_array.begin() ; it != index_array.end() ; ++it) {
        if(lm == *it) {
          break;
        }
        i++;
      }

      V = A[i] * (Q_pklm*(1-fR));
      V += (hR_conv_A[i] * (Q_pklm*fR*T0));
      V *= U[k];
      
      sigma += V;
    }
    
    // Fast Fourier Transform of the sum
    ComplexArray sigma_fft(nt);
    fft[p]->compute(sigma,sigma_fft);
    // The product (sigma * conj(E_i)) is made in the frequency domain
    // Then we go back to the time domain with an Inverse Fast Fourier Transform
    // Finally we multiply the result with i
    sigma_fft *= conj(E_i[p]);
    ifft[p]->compute(sigma_fft,N[p]);
    N[p] *= complex<double>(0.,sign_gamma);
  }

  return N;
}


ComplexArray Nonlinearity_withRaman::convolution_with_hR(ComplexArray& A)
{ 
  ComplexArray result_fft(nt);
  fft[0]->compute(A,result_fft);
  result_fft *= hR_fft;

  ComplexArray result(A);
  ifft[0]->compute(result_fft,result);
  result *= delta;
  
  return result;
}


void Nonlinearity_withRaman::compute_delayed_Raman_response(const double tau[2])
{
  ComplexArray hR(nt);

  double a = (1./(tau[0]*tau[0]) + 1./(tau[1]*tau[1])) * tau[0];
  for(unsigned int i = 0 ; i < nt ; i++) {
    double t = T0 * (i * delta);
    double b = exp(-t/tau[1]) * sin(t/tau[0]);
    hR[i] = complex<double>((a*b),0.);
  }
  
  fft[0]->compute(hR,hR_fft);
}
