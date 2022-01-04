#include "writeMatFile.hpp"

using namespace std;

typedef vector<ComplexArray> MultipleComplexArrays;


void writeMatFile(const MultipleComplexArrays phi, double cpu_time, const char* filename)
{
  unsigned int M = phi.size();
  unsigned int N = phi[0].getSize();
  unsigned int MN = M*N;
  
  double* real_part = new double[MN];
  double* imag_part = new double[MN];

  for(unsigned int p = 0; p < M ; p++) {
    unsigned int j = p*N;
    for(unsigned int i = 0 ; i < N ; i++) {
      real_part[j] = phi[p][i].real();
      imag_part[j] = phi[p][i].imag();
      j++;
    }
  }

  mat_complex_split_t result;
  result.Re = real_part;
  result.Im = imag_part;
  
  // Write
  mat_t *matfp = NULL; // matfp contains pointer to MAT file or NULL on failure
  matfp = Mat_CreateVer(filename, NULL, MAT_FT_MAT5); // or MAT_FT_MAT4 / MAT_FT_MAT73
  
  size_t dim2d[2] = {N,M};
  matvar_t *matvar = Mat_VarCreate("output",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dim2d,&result,MAT_F_COMPLEX);
  Mat_VarWrite(matfp,matvar,MAT_COMPRESSION_NONE);
  
  delete[] real_part;
  delete[] imag_part;
  
  size_t dim[2] = {1,1};
  matvar_t *matvar2 = Mat_VarCreate("cpu_time", MAT_C_DOUBLE,MAT_T_DOUBLE,2,dim,&cpu_time,0);
  Mat_VarWrite(matfp,matvar2,MAT_COMPRESSION_NONE);
  
  Mat_VarFree(matvar);
  Mat_VarFree(matvar2);
  Mat_Close(matfp);
}
