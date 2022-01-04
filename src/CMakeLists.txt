add_library(MMPropagation
  SHARED
  MultimodePropagation.cpp
  LawsonRK.cpp
  RungeKutta.cpp
  Nonlinearity.cpp
  ComplexArray.cpp
  Sparse3DMatrix.cpp
  DiscreteFourierTransform.cpp
  readMatFile.cpp
  writeMatFile.cpp)
  
install(TARGETS MMPropagation
             DESTINATION lib)
	
file(
    GLOB
    headers
    *.hpp
)

install(FILES ${headers}
  DESTINATION include/${CMAKE_PROJECT_NAME})
