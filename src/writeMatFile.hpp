#ifndef WRITE_MAT_FILE_HPP
#define WRITE_MAT_FILE_HPP

#include <matio.h>
#include "ComplexArray.hpp"

typedef vector<ComplexArray> MultipleComplexArrays;


void writeMatFile(const MultipleComplexArrays phi, double cpu_time, const char* filename);

#endif
