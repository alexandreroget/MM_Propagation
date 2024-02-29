#include "MultimodePropagation_GPU.cuh"


MultimodePropagationGPU::MultimodePropagationGPU(const struct SimulationParameters& in) : MultimodePropagation(in) {
  setGPUThreadConfiguration();
}


MultimodePropagationGPU::~MultimodePropagationGPU() {
  cufftDestroy(forward_plan_many);
  cufftDestroy(inverse_plan_many);
}


void MultimodePropagationGPU::setGPUThreadConfiguration() {
  cudaDeviceProp deviceProp;
  int deviceId = 0;
  cudaGetDeviceProperties(&deviceProp, deviceId);

  n_threads = deviceProp.maxThreadsPerBlock;
  
  n_blocks = nt / n_threads;
  if ((nt % n_threads) != 0) n_blocks++;
}


void MultimodePropagationGPU::setFFTPlansMany() {
 int batch = M;
 int rank = 1;

 int nRows = nt;
 int n[1] = {nRows};

 int idist = nRows;
 int odist = nRows;

 int inembed[] = {nRows};
 int onembed[] = {nRows};

 int istride = 1;
 int ostride = 1;

 cufftPlanMany(&forward_plan_many, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_Z2Z, batch);
 cufftPlanMany(&inverse_plan_many, rank, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2Z, batch);
}


void MultimodePropagationGPU::computeLawsonRK(const unsigned int nz) {
  initLawsonRK();

  ComplexArraysContainerGPU psi_GPU(nt, M);
  psi_GPU.getDataFromCPU(psi);

  ComplexArraysContainerGPU E1_GPU(nt, M);
  E1_GPU.getDataFromCPU(E1);
  
  ComplexArraysContainerGPU E_GPU(nt, M, RK.s);
  for(unsigned int i = 0 ; i < RK.s ; i++) {
    E_GPU.getDataFromCPU(E[i], i);
  }
  for(unsigned int i = 0 ; i < nz ; i++) {
    // Compute $\psi_{n+1}^{p}$ 
    applyRungeKuttaMethod(psi_GPU, E_GPU, delta_z);
    
    computeExponentialOperatorsMatrix(E_GPU, E1_GPU);
    z += delta_z;
  }
  
  psi_GPU.sendDataToCPU(psi);
}


void MultimodePropagationGPU::computeExponentialOperatorsMatrix(ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& E1_GPU) {
  // Compute $E_{n+1,i}^{p} = E_{n,i}^{p} E_{1,i}^p$
  for(unsigned int i = 0 ; i < RK.s ; i++) {
    multiplyManyComplexArrays(E_GPU, E1_GPU, i, 0); 
  }
}


void MultimodePropagationGPU::applyRungeKuttaMethod(ComplexArraysContainerGPU& psi_GPU, const ComplexArraysContainerGPU& E_GPU, double h) {
  ComplexArraysContainerGPU N(nt, make_cuDoubleComplex(0., 0.), M, RK.s);

  ComplexArraysContainerGPU result(nt, M);
  // $\psi_{n+1}^{p} = psi^{p}$
  copyManyComplexArrays(result, psi_GPU);

  for(unsigned int i = 0 ; i < RK.s ; i++) {
    ComplexArraysContainerGPU psi_i(nt, M);
    copyManyComplexArrays(psi_i, psi_GPU);
    
    // $\psi_{n,i}^{p} = psi_n^{p} + h \sum_{j=0}^{i-1} (a_{i,j} N_{n,j})$
    if(i != 0.) {
      for(unsigned int j = 0 ; j < i ; j++) {
        if(RK.a[i][j] != 0.) {
          ComplexArraysContainerGPU sigma(nt, M);
          copyManyComplexArrays(sigma, N, 0, j);
          multiplyManyComplexArraysByScalar(sigma, make_cuDoubleComplex(RK.a[i][j]*h, 0.));
          addManyComplexArrays(psi_i, sigma);
        }
      }
    }
    
    // Compute $N(z_n + c_i h, \psi_{n,i})$
    computeNonlinearity(N, E_GPU, psi_i, i);
    
    // $\psi_{n+1}^{p} = \psi_{n+1}^{p} + b_i h N(z_n + c_i h, \psi_{n,i})$
    if(RK.b[i] != 0.) {
      ComplexArraysContainerGPU sigma(nt, M);
      copyManyComplexArrays(sigma, N, 0, i);
      multiplyManyComplexArraysByScalar(sigma, make_cuDoubleComplex(RK.b[i]*h, 0.));
      
      addManyComplexArrays(result, sigma);
    }
  }

  copyManyComplexArrays(psi_GPU, result);
}


void MultimodePropagationGPU::computeManyFFT(ComplexArraysContainerGPU& u) {
  cufftExecZ2Z(forward_plan_many, u, u, CUFFT_FORWARD);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::computeManyIFFT(ComplexArraysContainerGPU& u) {
  cufftExecZ2Z(inverse_plan_many, u, u, CUFFT_INVERSE);
  cudaDeviceSynchronize();
  multiplyManyComplexArraysByScalar(u, make_cuDoubleComplex(1./nt, 0.));
}


ComplexArraysContainerCPU MultimodePropagationGPU::getResult() {
  ComplexArraysContainerGPU Phi_GPU(nt, M);
  Phi_GPU.getDataFromCPU(psi);

  computeManyFFT(Phi_GPU);
  
  ComplexArraysContainerCPU Lawson_operator(nt, M);
  for(unsigned int p = 0 ; p < M ; p++) {
    Lawson_operator[p] = exp(L[p]*z);
  }
  
  ComplexArraysContainerGPU Lawson_operator_GPU(nt, M);
  Lawson_operator_GPU.getDataFromCPU(Lawson_operator);
  
  multiplyManyComplexArrays(Phi_GPU, Lawson_operator_GPU);
  
  computeManyIFFT(Phi_GPU);
  
  multiplyManyComplexArraysByScalar(Phi_GPU, make_cuDoubleComplex(1./conversion_factor, 0.));

  ComplexArraysContainerCPU Phi(nt, M);
  Phi_GPU.sendDataToCPU(Phi);

  return Phi;
}

// --------------------------------------------------------------------------------------------- //

MultimodePropagationGPU_RamanOFF::MultimodePropagationGPU_RamanOFF(const struct SimulationParameters& in) : MultimodePropagationGPU(in) {
  setFFTPlans();
}


void MultimodePropagationGPU_RamanOFF::setFFTPlans() {
  setFFTPlansMany();
}


void MultimodePropagationGPU_RamanOFF::computeNonlinearity(ComplexArraysContainerGPU& N, const ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& psi_i, const unsigned int i) {
  // Step 1: Compute $U^{p} = E_{n,i}^{p} \psi_{n,i}^{p}$
  ComplexArraysContainerGPU U(nt, M);
  copyManyComplexArrays(U, psi_i);
  computeManyFFT(U);
  multiplyManyComplexArrays(U, E_GPU, 0, i);
  computeManyIFFT(U);
  
  ComplexArraysContainerGPU U_conj(nt, M);
  computeComplexConjugate(U_conj, U);
  
  // Step 2: Compute $V^{l,m} = U^{l} \overline{U^{m}}$
  unsigned int n_pairs = unique_pairs.size();
  ComplexArraysContainerGPU V(nt, n_pairs);
  {
    unsigned int lm = 0;
    for(const auto& pair : unique_pairs) {
      unsigned int l = pair % M;
      unsigned int m = pair / M;
    
      copySingleComplexArray(V, U, lm, l);
      multiplySingleComplexArray(V, U_conj, lm, m);
          
      lm++;
    }
  }
    
  // Step 3: Compute $\sigma^{p} = \sum_{k,l,m} {\tilde{Q^{p}_{k,l,m}} (U^{k} V^{l,m}))}$
  ComplexArraysContainerGPU sigma(nt, make_cuDoubleComplex(0., 0.), M);
  
  for(unsigned int p = 0 ; p < M ; p++) {
    unsigned int num_nzv = Q_tilde[p].getNumberOfNonZeroValues();
    
    for(unsigned int j = 0 ; j < num_nzv ; j++) {
      NonZeroValue nzv = Q_tilde[p].getNonZeroValue(j);
      unsigned int k = nzv.index[0];
      unsigned int l = nzv.index[1];
      unsigned int m = nzv.index[2];
      double Q_pklm = nzv.value;
      
      auto it = std::find(unique_pairs.begin(), unique_pairs.end(), (l + m*M));
      unsigned int lm = std::distance(unique_pairs.begin(), it);
        
      // Compute $W^{k,l,m} = \tilde{Q^{p}_{k,l,m}} U^{k} V^{l,m}$
      ComplexArraysContainerGPU W(nt, 1);
      
      copySingleComplexArray(W, U, 0, k);
      multiplySingleComplexArray(W, V, 0, lm);
      multiplyManyComplexArraysByScalar(W, make_cuDoubleComplex(Q_pklm, 0.));
      
      addSingleComplexArray(sigma, W, p, 0);
    }
  }
  
  // Step 4: Compute $N_{n,i}^{p} = \ii\gamma \overline{E_{n,i}^{p}} \sigma^{p}$
  computeManyFFT(sigma);

  ComplexArraysContainerGPU E_i_conj(nt, M);
  computeComplexConjugate(E_i_conj, E_GPU, 0, i);
  multiplyManyComplexArrays(sigma, E_i_conj);
  
  computeManyIFFT(sigma);

  multiplyManyComplexArraysByScalar(sigma, make_cuDoubleComplex(sign_gamma, 0.));
  
  copyManyComplexArrays(N, sigma, i, 0);
}

// --------------------------------------------------------------------------------------------- //

MultimodePropagationGPU_RamanON::MultimodePropagationGPU_RamanON(const struct SimulationParameters& in) : MultimodePropagationGPU(in), fR(in.raman_proportion), hR_fft(nt, 1) {
  setFFTPlans();

  ComplexArraysContainerCPU hR_CPU(nt, 1);
  for(unsigned int i = 0 ; i < nt ; i++) {
    hR_CPU[0][i] = std::complex<double>(in.raman_response[i], 0.);
  }
  
  hR_fft.getDataFromCPU(hR_CPU);
  computeSingleFFT(hR_fft);
}


void MultimodePropagationGPU_RamanON::setFFTPlans() {
 setFFTPlansMany();
 cufftPlan1d(&plan_1d, nt, CUFFT_Z2Z, 1);
}


void MultimodePropagationGPU_RamanON::computeNonlinearity(ComplexArraysContainerGPU& N, const ComplexArraysContainerGPU& E_GPU, const ComplexArraysContainerGPU& psi_i, const unsigned int i) {
  // Step 1: Compute $U^{p} = E_{n,i}^{p} \psi_{n,i}^{p}$
  ComplexArraysContainerGPU U(nt, M);
  copyManyComplexArrays(U, psi_i);
  computeManyFFT(U);
  multiplyManyComplexArrays(U, E_GPU, 0, i);
  computeManyIFFT(U);
  
  ComplexArraysContainerGPU U_conj(nt, M);
  computeComplexConjugate(U_conj, U);
  
  // Step 2: Compute $V^{l,m} = U^{l} \overline{U^{m}}$ and $V_Raman^{l,m} = h_R^{T_0} \star (U^{l} \overline{U^{m}})$
  unsigned int n_pairs = unique_pairs.size();
  ComplexArraysContainerGPU V(nt, n_pairs);
  ComplexArraysContainerGPU V_Raman(nt, n_pairs);
  {
    unsigned int lm = 0;
    for(const auto& pair : unique_pairs) {
      unsigned int l = pair % M;
      unsigned int m = pair / M;
    
      copySingleComplexArray(V, U, lm, l);
      multiplySingleComplexArray(V, U_conj, lm, m);
    
      copySingleComplexArray(V_Raman, V, lm, lm);
      convolveWithRamanResponse(V_Raman, lm);
    
      lm++;
    }
  }  
  
  // Step 3: Compute $\sigma^{p} = \sum_{k,l,m} {\tilde{Q^{p}_{k,l,m}} U^{k} ((1 - f_R) V^{l,m} + (f_R T_0) V_Raman^{l,m})}$
  ComplexArraysContainerGPU sigma(nt, make_cuDoubleComplex(0., 0.), M);
  for(unsigned int p = 0 ; p < M ; p++) {
    unsigned int num_nzv = Q_tilde[p].getNumberOfNonZeroValues();
    
    for(unsigned int j = 0 ; j < num_nzv ; j++) {
      NonZeroValue nzv = Q_tilde[p].getNonZeroValue(j);
      unsigned int k = nzv.index[0];
      unsigned int l = nzv.index[1];
      unsigned int m = nzv.index[2];
      double Q_pklm = nzv.value;
      
      auto it = std::find(unique_pairs.begin(), unique_pairs.end(), (l + m*M));
      unsigned int lm = std::distance(unique_pairs.begin(), it);
      
      // Compute $W^{k,l,m} = \tilde{Q^{p}_{k,l,m}} U^{k} ((1 - f_R) V^{l,m} + (f_R T_0) V_Raman^{l,m})$
      ComplexArraysContainerGPU W(nt, 1);

      copySingleComplexArray(W, U, 0, k);
      multiplySingleComplexArray(W, V, 0, lm);
      multiplyManyComplexArraysByScalar(W, make_cuDoubleComplex(Q_pklm*(1-fR), 0.));

      addSingleComplexArray(sigma, W, p, 0);
      
      copySingleComplexArray(W, U, 0, k);
      multiplySingleComplexArray(W, V_Raman, 0, lm);
      multiplyManyComplexArraysByScalar(W, make_cuDoubleComplex(Q_pklm*(fR*T0), 0.));
      
      addSingleComplexArray(sigma, W, p, 0);
    }
  }
  
  // Step 4: Compute $N_{n,i}^{p} = \ii\gamma \overline{E_{n,i}^{p}} \sigma^{p}$
  computeManyFFT(sigma);

  ComplexArraysContainerGPU E_i_conj(nt, M);
  computeComplexConjugate(E_i_conj, E_GPU, 0, i);
  multiplyManyComplexArrays(sigma, E_i_conj);
  
  computeManyIFFT(sigma);

  multiplyManyComplexArraysByScalar(sigma, make_cuDoubleComplex(sign_gamma, 0.));
  
  copyManyComplexArrays(N, sigma, i, 0);
}


void MultimodePropagationGPU_RamanON::convolveWithRamanResponse(ComplexArraysContainerGPU& V_Raman, const unsigned int index) {
  ComplexArraysContainerGPU result(nt, 1);
  copySingleComplexArray(result, V_Raman, 0, index);

  computeSingleFFT(result);
  multiplyManyComplexArrays(result, hR_fft);
  computeSingleIFFT(result);
  
  multiplyManyComplexArraysByScalar(result, make_cuDoubleComplex(delta_t, 0.));
  
  copySingleComplexArray(V_Raman, result, index, 0);
}


void MultimodePropagationGPU_RamanON::computeSingleFFT(ComplexArraysContainerGPU& u) {
  cufftExecZ2Z(plan_1d, u, u, CUFFT_FORWARD);
}


void MultimodePropagationGPU_RamanON::computeSingleIFFT(ComplexArraysContainerGPU& u) {
  cufftExecZ2Z(plan_1d, u, u, CUFFT_INVERSE);
  multiplyManyComplexArraysByScalar(u, make_cuDoubleComplex(1./nt, 0.));
}

// --------------------------------------------------------------------------------------------- //

void MultimodePropagationGPU::copyManyComplexArrays(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, M, 1);

  copyManyComplexArraysKernel<<<grid, block>>>(arraysContainerA, columnIndexA, arraysContainerB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::copySingleComplexArray(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int rowIndexA, const unsigned int rowIndexB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, 1, 1);

  copySingleComplexArrayKernel<<<grid, block>>>(arraysContainerA, rowIndexA, columnIndexA, arraysContainerB, rowIndexB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::addManyComplexArrays(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, M, 1);

  addManyComplexArraysKernel<<<grid, block>>>(arraysContainerA, columnIndexA, arraysContainerB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::addSingleComplexArray(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int rowIndexA, const unsigned int rowIndexB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, 1, 1);

  addSingleComplexArrayKernel<<<grid, block>>>(arraysContainerA, rowIndexA, columnIndexA, arraysContainerB, rowIndexB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::multiplyManyComplexArrays(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, M, 1);

  multiplyManyComplexArraysKernel<<<grid, block>>>(arraysContainerA, columnIndexA, arraysContainerB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::multiplySingleComplexArray(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int rowIndexA, const unsigned int rowIndexB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, 1, 1);

  multiplySingleComplexArrayKernel<<<grid, block>>>(arraysContainerA, rowIndexA, columnIndexA, arraysContainerB, rowIndexB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::multiplyManyComplexArraysByScalar(ComplexArraysContainerGPU& arraysContainerA, const cuDoubleComplex complexScalar, const unsigned int columnIndexA) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, M, 1);
  
  multiplyManyComplexArraysByScalarKernel<<<grid, block>>>(arraysContainerA, columnIndexA, complexScalar, nt, M);
  cudaDeviceSynchronize();
}


void MultimodePropagationGPU::computeComplexConjugate(ComplexArraysContainerGPU& arraysContainerA, const ComplexArraysContainerGPU& arraysContainerB, const unsigned int columnIndexA, const unsigned int columnIndexB) {
  dim3 block(n_threads, 1, 1);
  dim3 grid(n_blocks, M, 1);
  
  computeComplexConjugateKernel<<<grid, block>>>(arraysContainerA, columnIndexA, arraysContainerB, columnIndexB, nt, M);
  cudaDeviceSynchronize();
}

