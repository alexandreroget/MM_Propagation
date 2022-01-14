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
  
    double t_start = omp_get_wtime();
    mm_propagation.computeLawsonRK();
    double t_stop = omp_get_wtime();
    double cpu_time = (t_stop - t_start);
  
    MultipleComplexArrays phi_out(in.n_modes,ComplexArray(in.nt));
    phi_out = mm_propagation.getResult();
  
    writeMatFile(phi_out,0.,"output_MMPropagation.mat");
  }
  
  return 0;
}
