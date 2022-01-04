#ifndef READ_MAT_FILE_HPP
#define READ_MAT_FILE_HPP

#include <matio.h>
#include "MultimodePropagation.hpp"

typedef vector<double> doubleArray;
typedef vector<ComplexArray> MultipleComplexArrays;


void readMatFile(struct SimulationParameters& in, const char* input_file);
void readData(mat_t *mat, matvar_t *matvar, char* variable_name, const double*& data, unsigned int& data_size);
void readComplexData(mat_t *mat, matvar_t *matvar, char* variable_name, mat_complex_split_t*& complex_data, unsigned int& data_size);
void setDispersionCoefficients(const double* data, struct SimulationParameters& in, unsigned int data_size);
void setCouplingCoefficients(const double* data, struct SimulationParameters& in, unsigned int data_size);
void setInitialFields(mat_complex_split_t* complex_data, struct SimulationParameters& in, unsigned int data_size);

#endif
