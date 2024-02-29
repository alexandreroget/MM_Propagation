/*
* mm_propagation.cpp
*
* Compute the light propagation in a multimode fiber
*
* Usage : from MATLAB
*         >> outputs = mm_propagation(inputs)
*
* This is a C++ MEX-file for MATLAB.
* Copyright 2017 The MathWorks, Inc.
*
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <memory>

#include "MultimodePropagation_CPU.hpp"

class MexFunction : public matlab::mex::Function {

private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
  matlab::data::ArrayFactory factory;
  
  bool raman;
  
  bool save;
  std::string savename;
  unsigned int n_save;
  
public:
  /* Constructor for the class. */
  MexFunction()
  {
    matlabPtr = getEngine();
    
    raman = 0;
    n_save = 1;
  }

  void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) 
  {
    checkArguments(outputs, inputs);
    matlab::data::StructArray const matlabStructArray = inputs[0];
    checkStructureElements(matlabStructArray);
    outputs[0] = mm_propagation(matlabStructArray);
  }
  
  /* ------------------------------------------------------------------------------------------- */
  matlab::data::StructArray mm_propagation(matlab::data::StructArray const & matlabStructArray) 
  {
    double t_start = omp_get_wtime();
  
    struct SimulationParameters in;
    setInputParameters(in,matlabStructArray);
    
    MultimodePropagation* solver = NULL;
    in.raman_proportion == 0. ? solver = new MultimodePropagationCPU_RamanOFF(in) : solver = new MultimodePropagationCPU_RamanON(in);

    unsigned int nz = in.n_steps / n_save;
    double delta_Z = in.fiber_length / in.n_steps;
    double Z = 0.;
    
    for(unsigned int i = 0 ; i < n_save - 1 ; i++) {    
      solver->computeLawsonRK(nz);
    
      ComplexArraysContainer Phi(in.nt, in.n_modes);
      Phi = solver->getResult();
      Z  += nz * delta_Z;
      
      double t_stop = omp_get_wtime();
      double cpu_time = (t_stop - t_start);
      
      matlab::data::StructArray result = setResult(Phi, Z, cpu_time);
      
      std::string matlab_filename = savename + "_" + std::to_string(i+1);
      saveResult(result, matlab_filename);
    }
    
    nz += (in.n_steps % n_save);
    
    solver->computeLawsonRK(nz);
    
    ComplexArraysContainer Phi_out(in.nt, in.n_modes);
    Phi_out = solver->getResult();
    Z  += nz * delta_Z;
    
    double t_stop = omp_get_wtime();
    double cpu_time = (t_stop - t_start);
    
    matlab::data::StructArray result = setResult(Phi_out, Z, cpu_time);
    
    if(save) {
      std::string matlab_filename = savename;
      if(n_save > 1)  {
        matlab_filename += "_" + std::to_string(n_save);
      }
    
      saveResult(result, matlab_filename);
    }

    delete solver;

    return result;
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void setInputParameters(struct SimulationParameters & in, matlab::data::StructArray const & matlabStructArray) {
    auto fields = matlabStructArray.getFieldNames();
    std::vector<std::string> fieldNames(fields.begin(), fields.end());
    
    matlab::data::TypedArray<double> field1 = matlabStructArray[0]["n_modes"];
    in.n_modes = field1[0];
    
    matlab::data::TypedArray<double> field2 = matlabStructArray[0]["fiber_length"];
    in.fiber_length = field2[0];
  
    matlab::data::TypedArray<double> field3 = matlabStructArray[0]["dispersion_coefficients"];
    setDispersionCoefficients(in, field3);
        
    matlab::data::TypedArray<double> field4 = matlabStructArray[0]["coupling_coefficients"];
    setCouplingCoefficients(in, field4);
    
    if(matlabStructArray[0]["initial_fields"].getType() == matlab::data::ArrayType::DOUBLE) {
      matlab::data::TypedArray<double> field5 = matlabStructArray[0]["initial_fields"];
      setInitialFields(in, field5);
    }
    else {
      matlab::data::TypedArray<std::complex<double>> field5 = matlabStructArray[0]["initial_fields"];
      setInitialFields(in, field5);
    }
    
    matlab::data::TypedArray<double> field6 = matlabStructArray[0]["time_window"];
    in.time_window = field6[0];
    
    matlab::data::TypedArray<double> field7 = matlabStructArray[0]["pulse_width"];
    in.pulse_width = field7[0];
    
    matlab::data::TypedArray<double> field8 = matlabStructArray[0]["n_steps"];
    in.n_steps = field8[0];
    
    matlab::data::TypedArray<double> field9 = matlabStructArray[0]["method_order"];
    in.method_order = field9[0];
    
    matlab::data::TypedArray<double> field10 = matlabStructArray[0]["nonlinearity_const"];
    in.nonlinearity_const = field10[0];
    
    if(raman) {
      matlab::data::TypedArray<double> optional_field1 = matlabStructArray[0]["raman_proportion"];
      in.raman_proportion = optional_field1[0];
      if(in.raman_proportion != 0.) {
        matlab::data::TypedArray<double> optional_field2 = matlabStructArray[0]["raman_response"];
        setRamanParameters(in,optional_field2);
      }
    }
    else {
      in.raman_proportion = 0.;
    }
    
    if(save) {
      matlab::data::TypedArray<char16_t> optional_field3 = matlabStructArray[0]["savename"];
    
      savename.reserve(optional_field3.getNumberOfElements());
      for(auto it = optional_field3.begin(); it != optional_field3.end(); ++it) {
        savename += *it;
      }
      
      if(n_save != 1) {
        matlab::data::TypedArray<double> optional_field4 = matlabStructArray[0]["n_save"];
        n_save = optional_field4[0];
      }
    }
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void setDispersionCoefficients(struct SimulationParameters & in, matlab::data::TypedArray<double> field) {
    std::vector<long unsigned int> dim = field.getDimensions();
    unsigned int n_modes = dim[0];
    unsigned int L = dim[1];
    in.dispersion_coefficients = std::vector<doubleArray>(n_modes, doubleArray(L));
    
    unsigned int i = 0;
    unsigned int p = 0;
    for (auto elem : field) {
      in.dispersion_coefficients[p][i] = elem;
      p++;
      
      if(p == n_modes) {
        p = 0;
        i++;
      }
    }
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void setCouplingCoefficients(struct SimulationParameters & in, matlab::data::TypedArray<double> field) {
    std::vector<long unsigned int> dim = field.getDimensions();
    unsigned int M = dim[0];
    in.coupling_coefficients = std::vector<Sparse3DArray>(M);

    unsigned int p = 0;
    unsigned int k = 0;
    unsigned int l = 0;
    unsigned int m = 0;
    for (auto elem : field) {
      if(elem != 0.) {
        struct NonZeroValue nzv = {k, l, m, elem};
        in.coupling_coefficients[p].addNonZeroValue(nzv);
      }
      p++;
      
      if(p == M) {
        p = 0;
        k++;
        if(k == M) {
          k = 0;
          l++;
          if(l == M) {
            l = 0;
            m++;
          }
        }
      }
    }
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void setInitialFields(struct SimulationParameters & in, matlab::data::TypedArray<double> field) {
    std::vector<long unsigned int> dim = field.getDimensions();
    unsigned int M = dim[0];
    in.nt = dim[1];
    in.initial_fields = ComplexArraysContainer(in.nt, M);
    
    unsigned int i = 0;
    unsigned int p = 0;
    for (auto elem : field) {
      in.initial_fields[p][i] = std::complex<double>(elem,0.);
      p++;
      
      if(p == M) {
        p = 0;
        i++;
      }
    }
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void setInitialFields(struct SimulationParameters & in, matlab::data::TypedArray<std::complex<double>> field) {
    std::vector<long unsigned int> dim = field.getDimensions();
    unsigned int M = dim[0];
    in.nt = dim[1];
    in.initial_fields = ComplexArraysContainer(in.nt, M);
    
    unsigned int i = 0;
    unsigned int p = 0;
    for (auto elem : field) {
      in.initial_fields[p][i] = elem;
      p++;
      
      if(p == M) {
        p = 0;
        i++;
      }
    }
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void setRamanParameters(struct SimulationParameters & in, matlab::data::TypedArray<double> field) {
    in.raman_response = doubleArray(in.nt, 0.);

    unsigned int i = 0;
    for (auto elem : field) {
      in.raman_response[i] = elem;
      i++;
    }
  }
  
  /* ------------------------------------------------------------------------------------------- */
  matlab::data::StructArray setResult(const ComplexArraysContainer& Phi_out, const double Z, const double cpu_time) {
    unsigned int M = Phi_out.getNumberOfArrays();
    unsigned int N = Phi_out.getArraySize();
  
    std::vector<std::string> fieldNames = {"Z", "fields", "cpu_time"};
    matlab::data::StructArray result = factory.createStructArray({1, 1}, fieldNames);
    
    matlab::data::TypedArray<double> field1 = factory.createArray<double>({1, 1});
    field1[0] = Z;
    
    matlab::data::TypedArray<std::complex<double>> field2 = factory.createArray<std::complex<double>>({M, N});
    unsigned int i = 0;
    unsigned int p = 0;    
    for (auto elem : field2) {
      elem = Phi_out.getValue(p, i);
      p++;
      
      if(p == M) {
        p = 0;
        i++;
      }
    }
    
    matlab::data::TypedArray<double> field3 = factory.createArray<double>({1, 1});
    field3[0] = cpu_time;
    
    result[0]["Z"] = field1;
    result[0]["fields"] = field2;
    result[0]["cpu_time"] = field3;   
    
    return result; 
  }
  
  /* ------------------------------------------------------------------------------------------- */
  void saveResult(matlab::data::StructArray result, std::string matlab_filename) {
    matlabPtr->setVariable(u"output", result, matlab::engine::WorkspaceType::GLOBAL);

    std::vector<matlab::data::Array> args({
      factory.createCharArray(matlab_filename),
      factory.createCharArray("output")});
        
    const size_t numReturned = 0;
    matlabPtr->feval(u"save",numReturned,args);  
  }  
    
  /* This function makes sure that user has provided structure as input,
     and is not expecting more than one output in results. */
  void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) 
  {
    if(inputs.size() != 1) {
      displayError("One input required.");
    }
    if(outputs.size() > 1) {
      displayError("Too many outputs specified.");
    }
    if(inputs[0].getType() != matlab::data::ArrayType::STRUCT) {
      displayError("Input must be a structure.");
    }
  }
  
  /* ------------------- Make sure that the passed structure has valid data -------------------- */
  void checkStructureElements(matlab::data::StructArray const & matlabStructArray)
  { 
    size_t n_fields = matlabStructArray.getNumberOfFields();
    auto fields = matlabStructArray.getFieldNames();
    size_t total_num_of_elements = matlabStructArray.getNumberOfElements();
    
        
    if(total_num_of_elements != 1) {
      displayError("...");
    }

    if(n_fields < 10) {
      displayError(
        "Required fields:\n"
        "- n_modes: integer scalar\n"
        "- fiber_length: double scalar\n"
        "- dispersion_coefficients: double array (size: n_modes-by-L)\n"
        "- coupling_coefficients: double array (size: n_modes^4)\n"
        "- initial_fields: double or complex array (size: n_modes-by-N)\n"
        "- time_window: double scalar\n"
        "- pulse_width: double scalar\n"
        "- n_steps: integer scalar\n"
        "- method_order: integer scalar\n"
        "- nonlinearity_const: double scalar\n"
        "Optional fields:\n"
        "- raman_proportion: double scalar\n"
        "- raman_response: double array (size: N)"
        "- savename: char array\n"
        "- n_save: integer scalar\n");
    }
    
    std::vector<std::string> fieldNames(fields.begin(), fields.end());
    std::vector<std::string> requiredFields = {
      "n_modes", 
      "fiber_length", 
      "dispersion_coefficients",
      "coupling_coefficients",
      "initial_fields",
      "time_window",
      "pulse_width",
      "n_steps",
      "method_order",
      "nonlinearity_const"};

    for(unsigned int i = 0 ; i < 10 ; i++) {
      if(findFieldName(requiredFields[i], fieldNames, n_fields) == 0) {
        displayError("Field " + requiredFields[i] + " is missing.");
      }
    }
    
    if(matlabStructArray[0]["n_modes"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["n_modes"].getNumberOfElements() != 1) {
      invalidFieldInformation("n_modes");
      displayError("This field must contain an integer scalar.");
    }
    
    matlab::data::TypedArray<double> n_modes = matlabStructArray[0]["n_modes"];
    int M = n_modes[0];
    
    if(matlabStructArray[0]["fiber_length"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["fiber_length"].getNumberOfElements() != 1) {
      invalidFieldInformation("fiber_length");
      displayError("This field must contain a double scalar.");
    }
    
    if(matlabStructArray[0]["dispersion_coefficients"].getType() != matlab::data::ArrayType::DOUBLE) {
      invalidFieldInformation("dispersion_coefficients");
      displayError("This field must contain a n_modes-by-L double array.");
    }
    else {
      std::vector<long unsigned int> dim = matlabStructArray[0]["dispersion_coefficients"].getDimensions();
      if(dim.size() != 2 || dim[0] != M) {
        invalidFieldInformation("dispersion_coefficients");
        displayError("This field must contain a n_modes-by-L double array.");
      }
    }
    
    if(matlabStructArray[0]["coupling_coefficients"].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
      invalidFieldInformation("coupling_coefficients");
      displayError("This field must contain a n_modes^4 double array.");
    }
    else {
      if(M > 1) {
        std::vector<long unsigned int> dim = matlabStructArray[0]["coupling_coefficients"].getDimensions();
        if((dim.size() != 4) || ((dim[0] != M) || (dim[1] != M) || (dim[2] != M) || (dim[3] != M))) {
          invalidFieldInformation("coupling_coefficients");
          displayError("This field must contain a n_modes^4 double array.");
        }
      }
    }
    
    unsigned int N;

    if((matlabStructArray[0]["initial_fields"].getType() != matlab::data::ArrayType::DOUBLE) 
        && (matlabStructArray[0]["initial_fields"].getType() != matlab::data::ArrayType::COMPLEX_DOUBLE)) {
      invalidFieldInformation("initial_fields");
      displayError("This field must contain a n_modes-by-N double or complex array.");
    }
    else {
      std::vector<long unsigned int> dim = matlabStructArray[0]["initial_fields"].getDimensions();
      if(dim.size() != 2 || dim[0] != M) {
        invalidFieldInformation("initial_fields");
        displayError("This field must contain a n_modes-by-N double or complex array.");
      }
      
      N = dim[1];
    }
    
    if(matlabStructArray[0]["time_window"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["time_window"].getNumberOfElements() != 1) {
      invalidFieldInformation("time_window");
      displayError("This field must contain a double scalar.");
    }
    
    if(matlabStructArray[0]["pulse_width"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["pulse_width"].getNumberOfElements() != 1) {
      invalidFieldInformation("pulse_width");
      displayError("This field must contain a double scalar.");
    }
    
    if(matlabStructArray[0]["n_steps"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["n_steps"].getNumberOfElements() != 1) {
      invalidFieldInformation("n_steps");
      displayError("This field must contain an integer scalar.");
    }
    
    if(matlabStructArray[0]["method_order"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["method_order"].getNumberOfElements() != 1) {
      invalidFieldInformation("method_order");
      displayError("This field must contain an integer scalar.");
    }
    
    if(matlabStructArray[0]["nonlinearity_const"].getType() != matlab::data::ArrayType::DOUBLE
       || matlabStructArray[0]["nonlinearity_const"].getNumberOfElements() != 1) {
      invalidFieldInformation("nonlinearity_const");
      displayError("This field must contain a double scalar.");
    }
    
    // Check optional fields
    if(n_fields > 10) {    
      raman = findFieldName("raman_proportion", fieldNames, n_fields);
      
      if(raman) {
        if(findFieldName("raman_response", fieldNames, n_fields) == 0) {
          displayError("Field raman_response is missing.");
        }
        
        if(matlabStructArray[0]["raman_proportion"].getType() != matlab::data::ArrayType::DOUBLE 
           || matlabStructArray[0]["raman_proportion"].getNumberOfElements() != 1) {
          invalidFieldInformation("raman_proportion");
          displayError("This field must contain a double scalar.");
        }
    
        if((matlabStructArray[0]["raman_response"].getType() != matlab::data::ArrayType::DOUBLE) 
           && (matlabStructArray[0]["raman_response"].getType() != matlab::data::ArrayType::COMPLEX_DOUBLE)) {
            invalidFieldInformation("raman_response");
            displayError("This field must contain a N double array.");
        }
        else {
          std::vector<long unsigned int> dim = matlabStructArray[0]["raman_response"].getDimensions();
          if(dim[0] != N) {
            invalidFieldInformation("raman_response");
            displayError("This field must contain a N double array.");
          }
        }
      }

      save = findFieldName("savename", fieldNames, n_fields);
      
      if(save) {
        if(matlabStructArray[0]["savename"].getType() != matlab::data::ArrayType::CHAR) {
          invalidFieldInformation("savename");
          displayError("This field must contain a char array.");
        }
      
        if(findFieldName("n_save", fieldNames, n_fields)) {
          if(matlabStructArray[0]["n_save"].getType() != matlab::data::ArrayType::DOUBLE
             || matlabStructArray[0]["n_save"].getNumberOfElements() != 1) {
            invalidFieldInformation("n_save");
            displayError("This field must contain an integer scalar.");
          }
          
          n_save = 0;
        }
      }
    }
  }
  
  /* Find a specific field name in the passed structure */
  bool findFieldName(std::string field, std::vector<std::string> fields_array, unsigned int n_fields)
  {
    unsigned int i = 0;
    while((fields_array[i] != field) && (i < n_fields)) {
      i++;
    }
    return i < n_fields; 
  }
  
  /* Helper function to print output string on MATLAB command prompt. */
  void displayOnMATLAB(std::ostringstream stream)
  {
    matlabPtr->feval(u"fprintf", 0, std::vector<matlab::data::Array>({ factory.createScalar(stream.str())}));
  }
  
  /* Helper function to generate an error message from given string,
     and display it over MATLAB command prompt. */
  void displayError(std::string errorMessage)
  {
    matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({
      factory.createScalar(errorMessage) }));
  }
  
  /* Helper function to information about an empty field in the structure. */
  void invalidFieldInformation(std::string fieldName)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" is empty."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
  
  /* Helper function to information about an invalid field in the structure. */
  void emptyFieldInformation(std::string fieldName)
  {
    std::ostringstream stream;
    stream<<"Field: "<<std::string(fieldName)<<" contains wrong value."<<std::endl;
    displayOnMATLAB(std::move(stream));
  }
};
