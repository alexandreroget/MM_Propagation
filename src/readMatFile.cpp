#include "readMatFile.hpp"

using namespace std;

typedef vector<double> doubleArray;
typedef vector<ComplexArray> MultipleComplexArrays;


void readMatFile(struct SimulationParameters& in, const char* input_file)
{
  mat_t *mat = Mat_Open(input_file, MAT_ACC_RDONLY);
  
  if(mat) {
    matvar_t *matvar = NULL;
    const double *data = NULL;
    mat_complex_split_t* complex_data = NULL;
    unsigned int data_size;
    
    readData(mat,matvar,(char*) "n_modes",data,data_size);
    in.n_modes = (int) data[0];
    cout<<"n_modes = "<<in.n_modes<<endl;
  
    readData(mat,matvar,(char*) "fiber_length",data,data_size);
    in.fiber_length = data[0];
    cout<<"fiber_length = "<<in.fiber_length<<endl;
  
    readData(mat,matvar,(char*) "max_dispersion_order",data,data_size);
    in.max_dispersion_order = (int) data[0];
  
    readData(mat,matvar,(char*) "dispersion_coefficients",data,data_size);
    setDispersionCoefficients(data,in,data_size);

    readData(mat,matvar,(char*) "coupling_coefficients",data,data_size);
    setCouplingCoefficients(data,in,data_size);

    readData(mat,matvar,(char*) "N",data,data_size);
    in.nt = (int) data[0];

    readData(mat,matvar,(char*) "time_window",data,data_size);
    in.time_window = data[0];

    readData(mat,matvar,(char*) "pulse_width",data,data_size);
    in.pulse_width = data[0];

    readComplexData(mat,matvar,(char*) "initial_fields",complex_data,data_size);
    setInitialFields(complex_data,in,data_size);

    readData(mat,matvar,(char*) "n_steps",data,data_size);
    in.n_steps = data[0];

    readData(mat,matvar,(char*) "order_RK",data,data_size);
    in.order_RK = data[0];

    readData(mat,matvar,(char*) "nonlinearity_const",data,data_size);
    in.nonlinearity_const = data[0];
  
    readData(mat,matvar,(char*) "raman_proportion",data,data_size);
    in.raman_proportion = data[0];

    if(in.raman_proportion != 0.) {
      readData(mat,matvar,(char*) "raman_parameters",data,data_size);
      in.raman_parameters[0] = data[0];
      in.raman_parameters[1] = data[1];
    }
    else {
      in.raman_parameters[0] = 0.;
      in.raman_parameters[1] = 0.;
    }
    
    Mat_VarFree(matvar);
    Mat_Close(mat);
  }
  else {
    cerr<<"Error while opening the file "<<input_file<<"."<<endl;
    exit(EXIT_FAILURE);
  }
}


void readData(mat_t *mat, matvar_t *matvar, char* variable_name, const double*& data, unsigned int& data_size)
{
  matvar = Mat_VarRead(mat,variable_name);
  if(matvar) {
    data = static_cast<const double*>(matvar->data);
    data_size = matvar->nbytes/matvar->data_size;
  }
  else {
    cerr<<variable_name<<" not found."<<endl;
    Mat_VarFree(matvar);
    Mat_Close(mat);
    exit(EXIT_FAILURE);
  }
}


void readComplexData(mat_t *mat, matvar_t *matvar, char* variable_name, mat_complex_split_t*& complex_data, unsigned int& data_size)
{
  matvar = Mat_VarRead(mat,variable_name);
  if(matvar) {
    complex_data = (mat_complex_split_t*) matvar->data;
    data_size = matvar->nbytes/matvar->data_size;
  }
  else {
    cerr<<variable_name<<" not found."<<endl;
    Mat_VarFree(matvar);
    Mat_Close(mat);
    exit(EXIT_FAILURE);
  }
}


void setDispersionCoefficients(const double* data, struct SimulationParameters& in, unsigned int data_size)
{
  if(data_size == (in.n_modes * (in.max_dispersion_order+1))) {
    in.dispersion_coefficients = vector<doubleArray>(in.max_dispersion_order+1,doubleArray(in.n_modes));
  
    for(unsigned int i = 0 ; i <= in.max_dispersion_order ; i++) {
      for(unsigned int p = 0 ; p < in.n_modes ; p++) {
        unsigned int j = i*in.n_modes + p;
        in.dispersion_coefficients[i][p] = data[j];
      }
    }
  }
  else {
    cerr<<"The dispersion coefficients array size must be equal to "<<in.n_modes<<" x "<<(in.max_dispersion_order + 1)<<"."<<endl;
    exit(EXIT_FAILURE);
  }
}


void setCouplingCoefficients(const double* data, struct SimulationParameters& in, unsigned int data_size)
{
  if(data_size == pow(in.n_modes,4)) {
    in.coupling_coefficients = vector<Sparse3DMatrix>(in.n_modes);
  
    unsigned int M = in.n_modes;
    unsigned int M2 = M * M;
    unsigned int M3 = M2 * M;
  
    for(unsigned int p = 0 ; p < M ; p++) {
      for(unsigned int k = 0 ; k < M ; k++) {
        for(unsigned int l = 0 ; l < M ; l++) {
          for(unsigned int m = 0 ; m < M ; m++) {
            unsigned int i = m*M3 + l*M2 + k*M + p;
            if(data[i] != 0.) {
              in.coupling_coefficients[p].addNonZeroValue(k,l,m,data[i]);
            }
          }
        }
      }
    }
  }
  else {
    cerr<<"The coupling coefficients array size must be equal to "<<in.n_modes<<"^4."<<endl;
    exit(EXIT_FAILURE);
  }
}


void setInitialFields(mat_complex_split_t* complex_data, struct SimulationParameters& in, unsigned int data_size)
{
  if(data_size == (in.n_modes * in.nt)) {
    in.signal = MultipleComplexArrays(in.n_modes,ComplexArray(in.nt));
  
    double* Re = (double*)(complex_data->Re);
    double* Im = (double*)(complex_data->Im);

    for(unsigned int p = 0 ; p < in.n_modes ; p++) {
      for(unsigned int i = 0 ; i < in.nt ; i++) {
        unsigned int j = p*in.nt + i;
        in.signal[p][i] = complex<double>(Re[j],Im[j]);
      }
    }
  }
  else {
    cerr<<"The initial field array size must be equal to "<<in.n_modes<<" x "<<in.nt<<"."<<endl;
    exit(EXIT_FAILURE);
  }
}

