#include "MultimodePropagation.hpp"

typedef vector<double> doubleArray;
typedef vector<ComplexArray> MultipleComplexArrays;


MultimodePropagation::MultimodePropagation(struct SimulationParameters in) : 
M(in.n_modes), nt(in.nt), time_window(in.time_window), nz(in.n_steps), 
phi_out(in.n_modes,ComplexArray(in.nt,complex<double>(0.,0.)))
{
  double T0 = in.pulse_width;
  double t_final = time_window/T0;
  
  double gamma = in.nonlinearity_const;
    
  double max_beta2 = get_max(in.dispersion_coefficients[2]);
  double max_Q = get_max(in.coupling_coefficients);

  double LD = (T0*T0)/max_beta2;
  z_final = in.fiber_length/LD;

  conversion_factor = sqrt(abs(gamma)*LD) * sqrt(max_Q);

  vector<doubleArray> beta_tilde;
  set_beta_tilde(beta_tilde,in.dispersion_coefficients,T0,max_beta2);
 
  vector<Sparse3DMatrix> Q_tilde;
  set_Q_tilde(Q_tilde,in.coupling_coefficients,max_Q);
    
  for(unsigned int p = 0 ; p < M ; p++) {
    phi_in.push_back(in.signal[p]);
    phi_in[p] *= conversion_factor;
  }

  solver = new LawsonRK(M,
                        T0,
                        nt,
                        t_final,
                        in.order_RK,
                        gamma,
                        beta_tilde,
                        Q_tilde,
                        in.raman_proportion,
                        in.raman_parameters);
}


void MultimodePropagation::computeLawsonRK(unsigned int n)
{
  double delta_z = z_final/nz;
  unsigned int nz2 = nz/n;
  
  solver->initializeLawson(phi_in,delta_z);
  
  for(unsigned int i = 0 ; i < n ; i++) {
    MultipleComplexArrays phi_out(M,ComplexArray(nt));
    
    double z_stop = (i+1)*(z_final/n);
    phi_out = solver->compute(delta_z,nz2,z_stop);
    
    for(unsigned int p = 0 ; p < M ; p++) {
      phi_out[p] /= conversion_factor;
    }
    // saveResults(...);
  }
}


void MultimodePropagation::set_nz(unsigned int n_steps)
{
  nz = n_steps;
}


MultipleComplexArrays MultimodePropagation::getResult() const
{
  return phi_out;
}


void MultimodePropagation::saveResults(const MultipleComplexArrays phi, const string filename) const
{
  ofstream file;
  file.open(filename,ios::out|ios::trunc);

  double delta_t = (double) time_window/nt;
  
  double t = 0.;
  for(unsigned int i = 0 ; i < nt ; i++) {
    file<<t;
    
    for(unsigned int p = 0 ; p < M ; p++) {
      complex<double> v = phi[p][i];
      
      if(v.imag() >= 0) {
        file<<" ; "<<v.real()<<"+"<<v.imag()<<"j";
      }
      else {
        file<<" ; "<<v.real()<<v.imag()<<"j";
      }
    }
    
    t += delta_t;
    file<<endl;
  }
  
  file.close();
}


double MultimodePropagation::get_max(const vector<double> beta2) const
{
  double max_beta2 = 0.;
  
  for(unsigned int p = 0 ; p < M ; p++) {
    if(abs(beta2[p]) > max_beta2) {
      max_beta2 = abs(beta2[p]);
    }
  }
  return max_beta2;
}


double MultimodePropagation::get_max(const vector<Sparse3DMatrix> Q) const
{
  double max_Q = 0.;
  
  for(unsigned int p = 0 ; p < M ; p++) {
    unsigned int n_couplings = Q[p].getSize();
    for(unsigned int i = 0 ; i < n_couplings ; i++) {
      double Q_pklm = Q[p].getNonZeroValue(i);
      if(abs(Q_pklm) > max_Q) {
        max_Q = abs(Q_pklm);
      }
    }
  }
  return max_Q;
}


void MultimodePropagation::set_beta_tilde(vector<doubleArray>& beta_tilde, 
const vector<doubleArray> beta, const double T0, const double max_beta2)
{ 
  unsigned int n_beta = beta.size();
  double LD = (T0*T0)/max_beta2;
  
  doubleArray c;
  for(unsigned int i = 0 ; i < n_beta ; i++) {
    c.push_back(1./pow(T0,i) * LD);
  }

  for(unsigned int p = 0 ; p < M ; p++) {
    doubleArray beta_tilde_p;
    for(unsigned int i = 0 ; i < n_beta ; i++) {
      beta_tilde_p.push_back(beta[i][p]*c[i]); 
    }
    beta_tilde.push_back(beta_tilde_p);
  }
}


void MultimodePropagation::set_Q_tilde(vector<Sparse3DMatrix>& Q_tilde, const vector<Sparse3DMatrix> Q, const double max_Q)
{
  unsigned int k, l, m;
  double Q_pklm;

  for(unsigned int p = 0 ; p < M ; p++) {
    Sparse3DMatrix Q_tilde_p;
  
    unsigned int n_couplings = Q[p].getSize();
    for(unsigned int i = 0 ; i < n_couplings ; i++) {
      Q[p].getNonZeroValue(i,k,l,m,Q_pklm);
      Q_tilde_p.addNonZeroValue(k,l,m,Q_pklm/max_Q);
    }
    
    Q_tilde.push_back(Q_tilde_p);
  }
}
