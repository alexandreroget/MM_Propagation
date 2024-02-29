#include "MultimodePropagation_CPU.hpp"


MultimodePropagationCPU::MultimodePropagationCPU(const struct SimulationParameters& in) : MultimodePropagation(in) {
  setFFTPlans();
}


void MultimodePropagationCPU::setFFTPlans() {
  fft = new Forward_DFT*[M];
  ifft = new Inverse_DFT*[M];
    
  unsigned int max_threads = omp_get_max_threads();
  unsigned int n_threads = max_threads/M;
  if(n_threads == 0) n_threads++;
  
  for(unsigned int p = 0 ; p < M ; p++) {
    fft[p] = new Forward_DFT(n_threads,nt);
    ifft[p] = new Inverse_DFT(n_threads,nt);
  }
}


void MultimodePropagationCPU::computeLawsonRK(const unsigned int nz) {
  initLawsonRK();
  
  for(unsigned int i = 0 ; i < nz ; i++) {
    // Compute $\psi_{n+1}^{p}$
    applyRungeKuttaMethod(delta_z);
    
    computeExponentialOperatorsMatrix();
    z += delta_z;
  }
}


void MultimodePropagationCPU::computeExponentialOperatorsMatrix() {
  // Compute $E_{n+1,i}^{p} = E_{n,i}^{p} E_{1,i}^p$
  for(unsigned int i = 0 ; i < RK.s ; i++) {
    for(unsigned int p = 0 ; p < M ; p++) {
      E[i][p] *= E1[p];
    }
  }
}


void MultimodePropagationCPU::applyRungeKuttaMethod(const double h) {
  std::vector<ComplexArraysContainer> N(RK.s, ComplexArraysContainer(nt, std::complex<double>(0., 0.), M));
  
  // $\psi_{n+1}^{p} = psi^{p}$
  ComplexArraysContainer result(nt, M);
  for(unsigned int p = 0 ; p < M ; p++) {
    result[p] = psi[p];
  }
  
  for(unsigned int i = 0 ; i < RK.s ; i++) {
    // $\psi_{n,i}^{p} = psi_n^{p} + h \sum_{j=0}^{i-1} (a_{i,j} N_{n,j})$
    ComplexArraysContainer psi_i(nt, M);
    for(unsigned int p = 0 ; p < M ; p++) {
      psi_i[p] = psi[p];
      
      if(i != 0.) {
        for(unsigned int j = 0 ; j < i ; j++) {
          if(RK.a[i][j] != 0.) {
            psi_i[p] += (N[j][p]*(RK.a[i][j]*h));
          }
        }
      }
    }
    
    // Compute $N(z_n + c_i h, \psi_{n,i})$
    computeNonlinearity(N[i], E[i], psi_i);

    // $\psi_{n+1}^{p} = \psi_{n+1}^{p} + b_i h N(z_n + c_i h, \psi_{n,i})$
    if(RK.b[i] != 0.) {
      for(unsigned int p = 0 ; p < M ; p++) {
        result[p] += N[i][p]*(RK.b[i]*h);
      }
    }
  }
  
  for(unsigned int p = 0 ; p < M ; p++) {
    psi[p] = result[p];
  }
}


ComplexArraysContainer MultimodePropagationCPU::getResult() {
  ComplexArraysContainer Phi(nt, M);
  
  for(unsigned int p = 0 ; p < M ; p++) {
    ComplexArray psi_fft(nt);
    fft[0]->compute(psi[p], psi_fft);
      
    psi_fft *= exp(L[p]*z);
    
    ComplexArray phi(nt);
    ifft[0]->compute(psi_fft, phi);
      
    Phi[p] = phi / conversion_factor;
  }

  return Phi;
}

// --------------------------------------------------------------------------------------------- //

MultimodePropagationCPU_RamanOFF::MultimodePropagationCPU_RamanOFF(const struct SimulationParameters& in) : MultimodePropagationCPU(in) {}


void MultimodePropagationCPU_RamanOFF::computeNonlinearity(ComplexArraysContainer& N_i, ComplexArraysContainer& E_i, ComplexArraysContainer& psi_i) {
  // Step 1: Compute $U^{p} = E_{n,i}^{p} \psi_{n,i}^{p}$
  ComplexArraysContainer U(nt, M);
  #pragma omp parallel for shared(psi_i, E_i, U, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {
    ComplexArray U_fft(nt);
    fft[p]->compute(psi_i[p], U_fft);
    U_fft *= E_i[p];
    ifft[p]->compute(U_fft, U[p]);
  }
  
  // Step 2: Compute $V^{l,m} = U^{l} \overline{U^{m}}$
  unsigned int n_pairs = unique_pairs.size();
  ComplexArraysContainer V(nt, n_pairs);
  {
    unsigned int lm = 0;
    for(const auto& pair : unique_pairs) {
      unsigned int l = lm % M;
      unsigned int m = lm / M;

      V[lm] = U[l] * conj(U[m]);
      lm++;
    }
  }
  
  #pragma omp parallel for shared(Q_tilde, E_i, U, V, N_i, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {
    // Step 3: Compute $\sigma^{p} = \sum_{k,l,m} {\tilde{Q^{p}_{k,l,m}} (U^{k} V^{l,m}))}$
    ComplexArray sigma(nt, std::complex<double>(0., 0.));
    unsigned int num_nzv = Q_tilde[p].getNumberOfNonZeroValues();
    
    for(unsigned int j = 0 ; j < num_nzv ; j++) {
      NonZeroValue nzv = Q_tilde[p].getNonZeroValue(j);
      unsigned int k = nzv.index[0];
      unsigned int l = nzv.index[1];
      unsigned int m = nzv.index[2];
      double Q_pklm = nzv.value;
      
      auto it = std::find(unique_pairs.begin(), unique_pairs.end(), (l + m*M));
      unsigned int lm = std::distance(unique_pairs.begin(), it);
      
      ComplexArray W(nt);
      W = V[lm] * Q_pklm;
      W *= U[k];
         
      sigma += W;
    }
    
    // Step 4: Compute $N_{n,i}^{p} = \ii\gamma \overline{E_{n,i}^{p}} \sigma^{p}$
    ComplexArray sigma_fft(nt);
    fft[p]->compute(sigma, sigma_fft);
    sigma_fft *= conj(E_i[p]);
    ifft[p]->compute(sigma_fft, N_i[p]);
    N_i[p] *= std::complex<double>(0., sign_gamma);
  }
}

// --------------------------------------------------------------------------------------------- //

MultimodePropagationCPU_RamanON::MultimodePropagationCPU_RamanON(const struct SimulationParameters& in) : MultimodePropagationCPU(in), fR(in.raman_proportion), hR_fft(nt) {
  ComplexArray hR(nt);
  for(unsigned int i = 0 ; i < nt ; i++) {
    hR[i] = std::complex<double>(in.raman_response[i], 0.);
  }
  
  fft[0]->compute(hR, hR_fft);
}


void MultimodePropagationCPU_RamanON::computeNonlinearity(ComplexArraysContainer& N_i, ComplexArraysContainer& E_i, ComplexArraysContainer& psi_i) {
  // Step 1: Compute $U^{p} = E_{n,i}^{p} \psi_{n,i}^{p}$
  ComplexArraysContainer U(nt, M);
  #pragma omp parallel for shared(psi_i, E_i, U, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {
    ComplexArray U_fft(nt);
    fft[p]->compute(psi_i[p], U_fft);
    U_fft *= E_i[p];
    ifft[p]->compute(U_fft, U[p]);
  }
  
  // Step 2: Compute $V^{l,m} = U^{l} \overline{U^{m}}$ and $V_Raman^{l,m} = h_R^{T_0} \star (U^{l} \overline{U^{m}})$
  unsigned int n_pairs = unique_pairs.size();
  ComplexArraysContainer V(nt, n_pairs);
  ComplexArraysContainer V_Raman(nt, n_pairs);
  {
    unsigned int lm = 0;
    for(const auto& pair : unique_pairs) {
      unsigned int l = lm % M;
      unsigned int m = lm / M;

      V[lm] = U[l] * conj(U[m]);
      V_Raman[lm] = convolveWithRamanResponse(V[lm]);
      lm++;
    }
  }
  
  #pragma omp parallel for shared(Q_tilde, E_i, U, V, V_Raman, N_i, fft, ifft)
  for(unsigned int p = 0 ; p < M ; p++) {
    // Step 3: Compute $\sigma^{p} = \sum_{k,l,m} {\tilde{Q^{p}_{k,l,m}} U^{k} ((1 - f_R) V^{l,m} + (f_R T_0) V_Raman^{l,m})}$
    ComplexArray sigma(nt, std::complex<double>(0., 0.));
    unsigned int num_nzv = Q_tilde[p].getNumberOfNonZeroValues();
    
    for(unsigned int j = 0 ; j < num_nzv ; j++) {
      NonZeroValue nzv = Q_tilde[p].getNonZeroValue(j);
      unsigned int k = nzv.index[0];
      unsigned int l = nzv.index[1];
      unsigned int m = nzv.index[2];
      double Q_pklm = nzv.value;
      
      auto it = std::find(unique_pairs.begin(), unique_pairs.end(), (l + m*M));
      unsigned int lm = std::distance(unique_pairs.begin(), it);

      ComplexArray W(nt);
      W = V[lm] * (Q_pklm*(1-fR));
      W += (V_Raman[lm] * (Q_pklm*fR*T0));
      W *= U[k];
            
      sigma += W;
    }
    
    // Step 4: Compute $N_{n,i}^{p} = \ii\gamma \overline{E_{n,i}^{p}} \sigma^{p}$
    ComplexArray sigma_fft(nt);
    fft[p]->compute(sigma, sigma_fft);
    sigma_fft *= conj(E_i[p]);
    ifft[p]->compute(sigma_fft, N_i[p]);
    N_i[p] *= std::complex<double>(0., sign_gamma);
  }
}


ComplexArray MultimodePropagationCPU_RamanON::convolveWithRamanResponse(const ComplexArray& V) { 
  ComplexArray V_fft(nt);
  fft[0]->compute(V, V_fft);
  V_fft *= hR_fft;
    
  ComplexArray result(nt);
  ifft[0]->compute(V_fft, result);
  result *= delta_t;
    
  return result;
}

