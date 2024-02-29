#include "DiscreteFourierTransform.hpp"

DiscreteFourierTransform::DiscreteFourierTransform(unsigned int n_threads, unsigned int size, int direction) : n(size) {
  fftw_init_threads();

  omp_set_num_threads(n_threads);
  fftw_plan_with_nthreads(n_threads);
  
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

  plan = fftw_plan_dft_1d(n,in,out,direction,FFTW_ESTIMATE);
  
  fftw_free(in);
  fftw_free(out);
}
  

void DiscreteFourierTransform::destroy_plan() { fftw_destroy_plan(plan); }
  

Forward_DFT::Forward_DFT(unsigned int n_threads, unsigned int size) : DiscreteFourierTransform(n_threads, size, FFTW_FORWARD) {}


void Forward_DFT::compute(Complex* in, Complex* out) {
  fftw_execute_dft(plan, (fftw_complex*) in, (fftw_complex*) out);
}


Inverse_DFT::Inverse_DFT(unsigned int n_threads, unsigned int size) : DiscreteFourierTransform(n_threads, size, FFTW_BACKWARD) {}


void Inverse_DFT::compute(Complex* in, Complex* out) {
  fftw_execute_dft(plan, (fftw_complex*) in, (fftw_complex*) out);
  
  for(unsigned int i = 0 ; i < n ; i++) {
    out[i] *= std::complex<double>(1./n,0.);
  }
}
