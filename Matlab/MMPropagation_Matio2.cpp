#include "readMatFile.hpp"
#include "writeMatFile.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  char* input_filename = argv[1];

  if(input_filename == NULL) {
    cout<<"You must specify the .mat input filename."<<endl;
  }
  else {
    cout<<input_filename<<endl;
  
    struct SimulationParameters in;
    readMatFile(in,input_filename);    
    
    MultimodePropagation mm_propagation(in);
  
    MultipleComplexArrays phi_out(in.n_modes,ComplexArray(in.nt));

    for(unsigned int i = 0 ; i < 1 ; i++) {
      // unsigned int nz = (int) pow(2,(i+2));
      unsigned int nz = 8192;
      double t_start = omp_get_wtime();
      mm_propagation.computeLawsonRK(nz);
      double t_stop = omp_get_wtime();
      double cpu_time = (t_stop - t_start);

      phi_out = mm_propagation.getResult();
  
      char filename[100];
      sprintf(filename,"cpp_nz%d.mat",nz);
      writeMatFile(phi_out,cpu_time,filename);
    }
  }
  
  return 0;
}
