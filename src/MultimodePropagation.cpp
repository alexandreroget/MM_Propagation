#include "MultimodePropagation.hpp"


MultimodePropagation::MultimodePropagation(const struct SimulationParameters& in) : 
M(in.n_modes), nt(in.nt), T0(in.pulse_width), L(nt, M), psi(nt, M), E1(nt, M) {
  setRungeKuttaCoefficients(in.method_order);
  
  E = std::vector<ComplexArraysContainer>(RK.s, ComplexArraysContainer(nt, M));
  
  setDimensionlessVariables(in);
  setUniquePairs();
}


void MultimodePropagation::setRungeKuttaCoefficients(const unsigned int method_order) {
  switch(method_order) {
  case 8:
    Butcher11(RK);
    break;
  case 6:
    Butcher7(RK, 1.);
    break;
  default:
    Butcher4(RK);
    break;
  }
}


/* ----- Change of variables that allows to reduce to dimensionless variables and unknowns ----- */
void MultimodePropagation::setDimensionlessVariables(const struct SimulationParameters& in) {
  double gamma = in.nonlinearity_const;
  gamma >= 0 ? sign_gamma = 1 : sign_gamma = -1;
  
  std::vector<doubleArray> beta = in.dispersion_coefficients;    
  double max_beta2 = getMax(beta);

  std::vector<Sparse3DArray> Q = in.coupling_coefficients;
  double max_Q = getMax(Q);
  
  // Caracteristic length
  double LD = (T0*T0)/max_beta2;
  
  beta_tilde = setDimensionlessDispersionCoefficients(beta, LD, max_beta2);
  Q_tilde = setDimensionlessCouplingCoefficients(Q, max_Q);
    
  t_final = in.time_window/T0;
  delta_t = (double) t_final/nt;
  
  z_final = in.fiber_length/LD;
  delta_z = z_final/in.n_steps;
  
  conversion_factor = std::sqrt(std::abs(gamma)*LD) * std::sqrt(max_Q);

  setDimensionlessInitialFields(in.initial_fields);
}


void MultimodePropagation::setDimensionlessInitialFields(ComplexArraysContainer initial_fields) {
  for(unsigned int p = 0 ; p < M ; p++) {
    psi[p] = initial_fields[p] * conversion_factor;
  }
}


/* ----- Get $ max_p{|\beta_2^{(p)}|} $ ----- */
double MultimodePropagation::getMax(const std::vector<doubleArray>& beta) const {
  double max_beta2 = 0.;
  
  for(unsigned int p = 0 ; p < M ; p++) {
    if(std::abs(beta[p][2]) > max_beta2) {
      max_beta2 = std::abs(beta[p][2]);
    }
  }

  return max_beta2;
}
   

/* ----- Get $ \max_{k,l,m,p}{|Q_{k,l,m}^{(p)}|} $ ----- */
double MultimodePropagation::getMax(const std::vector<Sparse3DArray>& Q) const {
  double max_Q = 0.;
  
  for(unsigned int p = 0 ; p < M ; p++) {
    unsigned int num_nzv = Q[p].getNumberOfNonZeroValues();
    for(unsigned int i = 0 ; i < num_nzv ; i++) {
      struct NonZeroValue nzv = Q[p].getNonZeroValue(i);  
      double Q_pklm = nzv.value;
        
      if(std::abs(Q_pklm) > max_Q) {
        max_Q = std::abs(Q_pklm);
      }
    }
  }

  return max_Q;
}


std::vector<doubleArray> MultimodePropagation::setDimensionlessDispersionCoefficients(const std::vector<doubleArray>& beta, const double LD, const double max_beta2) {
  // Compute $\tilde{\beta_l^{(p)}} = (L_D / T_0^l) \beta_l^{(p)}$
  unsigned int n_beta = beta[0].size();
  std::vector<doubleArray> beta_tilde(M, doubleArray(n_beta));
      
  doubleArray c;
  for(unsigned int i = 0 ; i < n_beta ; i++) {
    c.push_back(1./pow(T0,i) * LD);
  }
  
  for(unsigned int p = 0 ; p < M ; p++) {
    for(unsigned int i = 0 ; i < n_beta ; i++) {
      beta_tilde[p][i] =  beta[p][i] * c[i];
    }
  }
  
  return beta_tilde;
}


std::vector<Sparse3DArray> MultimodePropagation::setDimensionlessCouplingCoefficients(const std::vector<Sparse3DArray>& Q, const double max_Q) {
  // Compute $\tilde{Q_{k,l,m}^{(p)}} = Q_{k,l,m}^{(p)} / \max_{k,l,m,p}{|Q_{k,l,m}^{(p)}|}$
  std::vector<Sparse3DArray> Q_tilde(M);
  
  for(unsigned int p = 0 ; p < M ; p++) { 
    unsigned int num_nzv = Q[p].getNumberOfNonZeroValues();
    for(unsigned int i = 0 ; i < num_nzv ; i++) {
      NonZeroValue nzv = Q[p].getNonZeroValue(i);
        
      nzv.value /= max_Q;
        
      Q_tilde[p].addNonZeroValue(nzv);
    }
  }
  
  return Q_tilde;
}


void MultimodePropagation::setUniquePairs() {
  for(unsigned int p = 0 ; p < M ; p++) {
    unsigned int num_nzv = Q_tilde[p].getNumberOfNonZeroValues();
    for(unsigned int i = 0 ; i < num_nzv ; i++) {
      NonZeroValue nzv = Q_tilde[p].getNonZeroValue(i);
        
      unsigned int l = nzv.index[1];
      unsigned int m = nzv.index[2];
      unsigned int lm = l + m*M;
      
      auto it = std::find(unique_pairs.begin(), unique_pairs.end(), lm);
      if(it == unique_pairs.end()) {
        unique_pairs.push_back(lm);
      }
    }
  }
}


void MultimodePropagation::initLawsonRK() {
  z = 0.;
  
  computeLinearOperators();
  initExponentialOperatorsMatrix();
}

void MultimodePropagation::computeLinearOperators() {
  // Compute $L_p = \ii L_D \sum_{\ell=0}^L \frac{\beta^{(p)}_\ell}{\ell !} (\frac{\ii}{T_0}\frac{\partial}{\partial t})^\ell$  
  std::complex<double> i = std::complex<double>(0., 1.);
  unsigned int n_beta = beta_tilde[0].size();
    
  std::vector<double> alpha(nt);
  for(unsigned int j = 0 ; j < nt ; j++) {
    int k;
    j<=nt/2 ? k=j : k=j-nt;
    alpha[j] = (2*M_PI/t_final) * k;
  }
  
  for(unsigned int p = 0 ; p < M ; p++) {
    for(unsigned int j = 0 ; j < nt ; j++) {
      L[p][j] = i * (beta_tilde[p][0] -beta_tilde[p][1]*alpha[j] -beta_tilde[p][2]/2*(alpha[j]*alpha[j]));
      if(n_beta >= 3) {
        for(unsigned int k = 3 ; k < n_beta ; k++) {
          L[p][j] += (i * (pow(-1,k) * beta_tilde[p][k]/factorial(k) * pow(alpha[j],k)));
        }
      }
    }
  }
}


void MultimodePropagation::initExponentialOperatorsMatrix() {
  for(unsigned int p = 0 ; p < M ; p++) {
    // Compute $E_{1,i}^{p} = {\rm e}^{h \mathbb L_p}$
    E1[p] = exp(L[p]*delta_z);
    // Compute $E_{0,i}^{p} = {\rm e}^{c_i h \mathbb L_p}$
    for(unsigned int i = 0 ; i < RK.s ; i++) {
      E[i][p] = exp(L[p]*(RK.c[i]*delta_z));
    }
  }
}


int MultimodePropagation::factorial(const int a) {
  int result = 1;
  if(a > 0) {
    for(int i = 1 ; i <= a ; i++) {
      result *= i;
    }
  }
  return result;
}


void MultimodePropagation::computeLawsonRK() {
  unsigned int nz = floor(z_final / delta_z);
  computeLawsonRK(nz);
}

