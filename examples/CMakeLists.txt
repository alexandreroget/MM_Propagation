include_directories(${CMAKE_SOURCE_DIR}/src)

file(
    GLOB
    usage_examples
    *.cpp
)

foreach(f ${usage_examples})
  get_filename_component(exampleName ${f} NAME_WE)
  add_executable(${exampleName} ${f})
  target_link_libraries(${exampleName} lib_MMPropagation)
  target_link_libraries(${exampleName} ${FFTW_LIBRARIES})
  target_link_libraries(${exampleName} ${MATIO_LIBRARIES})
  #target_link_libraries(${exampleName} PRIVATE OpenMP::OpenMP_CXX)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${exampleName}
  DESTINATION bin
  RENAME ${CMAKE_PROJECT_NAME}-${exampleName})
endforeach(f)
