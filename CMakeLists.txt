cmake_minimum_required(VERSION 2.8.2)

project(MMPropagation)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

find_package(FFTW REQUIRED)
find_package(OpenMP REQUIRED)

set(SOURCES
    src/MultimodePropagation_CPU.cpp
    src/MultimodePropagation.cpp
    src/ComplexArray.cpp
    src/ComplexArraysContainer.cpp
    src/Sparse3DArray.cpp
    src/DiscreteFourierTransform.cpp
    src/RungeKutta.cpp
)

include_directories(include)

add_library(mmpropagation STATIC ${SOURCES})

file(GLOB SIMULATIONS_SOURCES "simulations/*.cpp")

foreach(simulation_source ${SIMULATIONS_SOURCES})
    get_filename_component(simulation_name ${simulation_source} NAME_WE)
    
    add_executable(${simulation_name} ${simulation_source})
    
    target_link_libraries(${simulation_name} PRIVATE mmpropagation ${FFTW_LIBRARIES} OpenMP::OpenMP_CXX)
    target_compile_options(${simulation_name} PRIVATE "${OpenMP_CXX_FLAGS}")
    target_include_directories(${simulation_name} PUBLIC ${MMPropagation_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})
endforeach()
