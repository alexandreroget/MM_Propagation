include_directories(${CMAKE_SOURCE_DIR}/src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "../findFFTW")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "../findMATIO")
  
find_package(FFTW REQUIRED)
find_package(MATIO REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(mainMatio mainMatio.cpp)

target_link_libraries(mainMatio PRIVATE MMPropagation ${FFTW_LIBRARIES} ${MATIO_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_options(mainMatio PRIVATE "${OpenMP_CXX_FLAGS}")
target_include_directories(mainMatio PUBLIC ${FFTW_INCLUDE_DIRS})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mainMatio
    DESTINATION bin
    RENAME ${CMAKE_PROJECT_NAME})
