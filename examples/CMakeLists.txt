include_directories(${CMAKE_SOURCE_DIR}/src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "../findFFTW")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "../findMATIO")

find_package(FFTW REQUIRED)
find_package(MATIO REQUIRED)
find_package(OpenMP REQUIRED)

file(
    GLOB
    usage_examples
    *.cpp
)

foreach(f ${usage_examples})
  get_filename_component(exampleName ${f} NAME_WE)
  add_executable(${exampleName} ${f})
  
  target_link_libraries(${exampleName} PRIVATE MMPropagation ${FFTW_LIBRARIES} ${MATIO_LIBRARIES} OpenMP::OpenMP_CXX)
  target_compile_options(${exampleName} PRIVATE "${OpenMP_CXX_FLAGS}")
  target_include_directories(${exampleName} PUBLIC ${FFTW_INCLUDE_DIRS})

#  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${exampleName}
#  DESTINATION bin
#  RENAME ${CMAKE_PROJECT_NAME}-${exampleName})
endforeach(f)
