# OVERVIEW

This user guide describes how to use the C++ software MM_Propagation, developped in the Paradyse team of Inria by Alexandre Roget.

MM_Propagation is a C++ code that enables to simulation the light propagation through a multimode optical fiber.

Numerically, this amounts to solving a system of coupled nonlinear nonlocal Schrödinger evolution equations along the abscissa of the fiber.
In order achieve high efficiency (meaning required CPU time to achieve a given numerical error), we propose to use high order explicit Runge-Kutta methods, coupled with a Lawson transform.
This allows for fast numerical computations (we take advantage of Fast Fourier transforms and we use explicit methods that do not require CFL stability conditions, thanks to Lawson transformation), which achieve high precision (since the underlying Runge-Kutta methods have high order).

# REQUIRED LIBRARIES

## OpenMP

The computation of MM_Propagation requires the application programming interface **OpenMP** that supports multi-platform shared-memory multiprocessing. The parallel computing enables to gradually reduce the computational time.

If OpenMP is no featured in the compiler, you can configure it by running the following command in your terminal:

```sh
sudo apt-get install libomp-dev
```

## FFTW

The numerical one step method is coupled with a uniform time discretization and the action of dispersive operators
along the fiber is computed via fast Fourier transforms.

To do so, MM_Propagation also requires the **FFTW**, which is a C subroutine library for computing the Discrete Fourier Transform (DFT) in one or more dimensions, of arbitrary input size, and of both real and complex data. Get the latest stable release of FFTW on the download page: https://www.fftw.org/download.html

See also https://www.fftw.org/fftw3_doc/Installation-on-Unix.html

The **--enable-openmp** flag is necessary in order to induce parallelism. The installation commands are thus the following:

```sh
./configure --enable-openmp
make
sudo make install
```

### Installation on macOS

On macOS, you may encounter the following error message:

```sh
configure: error: C compiler cannot create executables
```

If this happens, please run the following command instead:

```sh
./configure CC=gfortran --enable-openmp CFLAGS="-fopenmp"
```

You can install gcc (which includes gfortran) with the following command:

```sh
brew install gcc
```

The installation is now complete.

# USING MM_PROPAGATION WITH A C++ CODE

First, you have to create the source file (.cpp) in the **simulations** folder. 

Simulation settings for a linear propagation in a GRIN fiber is also available in this folder. You can use this example as a reference for your own simulation.

## Input parameters

In order to run a simulation, you first have to create an object of `MM_Propagation`.
All required parameters are contained in a structure called `InputParameters`. Its fields are the following:

|Field                     |Name                      |Unit          |Data type                        |Size          | 
|:-------------------------|:-------------------------|:-------------|:--------------------------------|:-------------|
|Number of modes ($M$)     |n_modes                   |              |integer                          |$1$           |
|Fiber length              |fiber_length              |$m$           |double                           |$1$           |
|Dispersion coefficients   |dispersion_coefficients   |$ps^i/m$      |std::vector<std::vector<double>> |$M \times L$  |
|Coupling coefficients     |coupling_coefficients     |$m^{-2}$      |std::vector<Sparse3DArray>       |$M$           |
|Initial fields            |initial_fields            |$W^{\frac12}$ |`ComplexArraysContainer`         |$M \times NT$ |
|Time window               |time_window               |$ps$          |double                           |$1$           |
|Pulse width               |pulse_width               |$ps$          |double                           |$1$           |
|Number of steps           |n_steps                   |              |integer                          |$1$           |
|Method order              |method_order              |              |integer                          |$1$           |
|Nonlinearity constant     |nonlinearity_const        |$W^{-1}m$     |double                           |$1$           |
|Raman proportion          |raman_proportion          |              |double                           |$1$           |
|Raman response            |raman_response            |$ps$          |std::vector<double>              |$NT$          |

with $L$ the maximal dispersion order.

## Secondary classes

1. Sparse 3D Array

In order to save computational time and memory requirements, only non-zero coupling coefficients are saved. Thus the `Sparse3DArray` class has been created. This class enables to set non zero coefficients and their corresponding index $k$, $l$ and $m$. All these parameters are stored in a structure called `NonZeroValue`.

Here's an example:

```c++
double SR = 6.4341e12;
  
Sparse3DArray Q;
struct NonZeroValue nzv;
  
// Q_{k = 0, l = 0, m = 0}
nzv = {0,0,0,SR};
Q.addNonZeroValue(nzv);

// Q_{k = 0, l = 0, m = 1}
nzv = {0,0,1,SR/3};
Q.addNonZeroValue(nzv);

// Q_{k = 0, l = 1, m = 1}
nzv = {0,1,1,SR/2};
Q.addNonZeroValue(nzv);

// Q_{k = 1, l = 1, m = 1}
nzv = {1,1,1,-SR/3};
Q.addNonZeroValue(nzv);
```

`coupling_coefficients` is defined in a `std::vector<Sparse3DArray>` of size $M$. Each element of the vector contains the non-zero coupling coefficients of a mode.
For example, `coupling_coefficients[0]` is the sparse 3D array of $Q^{(0)}$.

2. Complex Arrays

The solver makes many computations involving complex arrays. Therefore, two classes have been created: `ComplexArray` and `ComplexArraysContainer`.
These classes use the `complex` class of the C++ standard library.
The operator overloading enables to make computations with complex arrays. It also enables to easily set and get values from an array.

The following example shows how to build two solitons:

```c++
unsigned int M = 2;
unsigned int nt = 1024;
ComplexArraysContainer Phi(nt, M);

std::complex<double> i(0., 1.);

double c = 0.5;
double q = 8.;
double t_final = 125.;

double h = (double) t_final/nt;
double t = -t_final/2;
double t0 = t_final/10;

for(unsigned int j = 0 ; j < nt ; j++) {
  phi[0][j] = 0.5 * sqrt(q/2) * (1./cosh(q/4 * (t - t0))) *  exp(i * c * ((t - t0)/2));
  phi[1][j] = 0.5 * sqrt(q/2) * (1./cosh(q/4 * (t + t0))) *  exp(i * c * ((t + t0)/2));
  t += h;
}
```

`initial_fields` is defined in a `ComplexArraysContainer` of size $M \times NT$. Each element of the vector contains the initial field of a mode.
For example, `initial_fields[0]` contains all values of $\Phi^{(0)}(Z=0,:)$.

## Run the propagation on CPU

If you wish to run the propagation without the Raman effect, use the `MultimodePropagationCPU_RamanOFF` class. Otherwise, use the `MultimodePropagationCPU_RamanON` class.

```c++
/* --- Set input parameters --- */
struct InputParameters input;
input.n_modes = 2;
input.fiber_length = 100;
// ...

/* ----- Create the class ----- */
MultimodePropagationCPU_RamanOFF solver(input);

/* ---- Run the propagation --- */
solver.computeLawsonRK();

/* ------ Get the result ------ */
ComplexArraysContainer output = solver.getResult();
```

### Compilation

The makefile is generated by CMake. Download the latest release on the following page: https://cmake.org/install/

Then, go to the main folder and write in a terminal the following commands:

```sh
mkdir build
cd build
cmake..
make 
```

In addition to compiling the entire MM_Propagation code, the makefile compiles all the source files in the **simulations** folder and generates the corresponding executable files. Once the compilation is complete, you can then run the simulation.

## Run the propagation on GPU

If you wish to run the propagation without the Raman effect, use the `MultimodePropagationGPU\_RamanOFF` class. Otherwise, use the `MultimodePropagationGPU\_RamanON` class.

```c++
/* --- Set input parameters --- */
struct InputParameters input;
input.n_modes = 2;
input.fiber_length = 100;
// ...

/* ----- Create the class ----- */
MultimodePropagationGPU_RamanOFF solver(input);

/* ---- Run the propagation --- */
solver.computeLawsonRK();

/* ------ Get the result ------ */
ComplexArraysContainer output = solver.getResult();
```

### Compilation

1. Put your CUDA code (.cu file) in the `simulations` directory.
2. Go to the main folder and execute the following command to compile the CUDA code:

```sh
nvcc -o MM_Propagation_GPU simulations/main_GPU.cu src/MultimodePropagation_GPU.cu src/MultimodePropagation.cpp src/ComplexArraysContainerGPU.cu src/ComplexArraysContainer.cpp src/ComplexArray.cpp src/RungeKutta.cpp src/Sparse3DArray.cpp -I include -std=c++11 -lcufft
```

Please ensure that **nvcc** (CUDA compiler) is installed on your system. 

The resulting executable is able to run the propagation simulation on the GPU.

# USING MM_PROPAGATION WITH MATLAB

You can also use MM_Propagation from MATLAB. To do so, go to the **matlab** folder.

## Creating the MEX File

A MEX file is a dynamically linked subroutine that is compiled from C or C++ source code and can be called from within MATLAB, just like any other MATLAB function. You can generate the MEX file with the following command:

```sh
mex CXXFLAGS="\$CXXFLAGS -std=c++11" CXXOPTIMFLAGS='\$CXXOPTIMFLAGS -Ofast -DNDEBUG' LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O3" -lgomp -lfftw3_omp -lfftw3 -lm -I../include mm_propagation.cpp ../src/*.cpp
```

### MEX File on MacOS

On macOS, no C compiler is supplied with MATLAB. To install MEX, you will need the Apple's development environment for macOS (Xcode) which is available on the Mac App Store.

Once Xcode is installed, run the following command in MATLAB:

```sh
mex -setup
```

You may see the following message:

Warning: Xcode is installed, but its license has not been accepted. Run Xcode and accept its license agreement.

In this case, you need to accept the Xcode license. To do so, run the following command:

```sh
sudo xcode-select --switch /Applications/Xcode.app/Contents/Developer
```

and then the following command:

```sh
sudo xcodebuild -license accept
```

## Run a simulation

Once the MEX file has been generated, you can use `mm_propagation` as a Matlab function. The syntax is the following:

```matlab
output = mm_propagation(input);
```

with `input` a structure containing ten required fields and four optional fields:

|  |Field                     |Name                     |Unit         |Data type                 |Size           | 
|--|:-------------------------|:------------------------|:------------|:-------------------------|:--------------|
|1 |Number of modes (M)       |n_modes                  |             |integer scalar            |$1$            |
|2 |Fiber length              |fiber_length             |$m$          |double scalar             |$1$            |
|3 |Dispersion coefficients   |dispersion_coefficients  |$ps^i/m$     |double array              |$M \times L$   |
|4 |Coupling coefficients     |coupling_coefficients    |$m^{-2}$     |double array              |$M^4$          |
|5 |Initial fields            |initial_fields           |$W^{\frac12}$|double or complex array   |$M \times NT$  |
|6 |Time window               |time_window              |$ps$         |double scalar             |$1$            |
|7 |Pulse width               |pulse_width              |$ps$         |double scalar             |$1$            |
|8 |Number of steps           |n_steps                  |             |integer scalar            |$1$            |
|9 |Method order              |method_order             |             |integer scalar            |$1$            |
|10|Nonlinearity constant     |nonlinearity_const       |$W^{-1}m$    |double scalar             |$1$            |
|11|Raman proportion          |raman_proportion         |             |double scalar             |$1$            |
|12|Raman response            |raman_response           |$ps$         |double array              |$NT$           |
|13|Save file name            |savename                 |             |char array                |$1$            |
|14|Number of saved files     |n_save                   |             |double scalar             |$1$            |

Fields 1 to 10 are mandatory. The field order doesn't matter, but the name of each field has to be respected. Otherwise, the function will not be able to read all input parameters correctly. Arrays sizes must be valid too. For example, if the `coupling_coefficients` array size is not equal to $M^4$, then an error will occur and the program will not compute the pulse propagation.

Fields 11 and 12 are optional. If `raman_proportion` is not specified, then the Raman scattering effect will not be taken into account during the computation. If `raman_proportion` is specified but not `raman_response`, then an error will occur.

Fields 13 and 14 are also optional. If you want to save results in a .mat file, then you have to specify `savename`.
The saved file (`savename.mat`) contains a structure of three fields:

|  |Field                 |Name     |Unit          |Data type     |Size           | 
|--|:---------------------|:--------|:-------------|:-------------|:--------------|
|1 |Position in the fiber |Z        |$m$           |double scalar |$1$            |
|2 |Fields                |fields   |$W^{\frac12}$ |complex array |$1$            |
|3 |CPU Time              |cpu_time |$s$           |double scalar |$M \times L$   |

You can specify the number of .mat files you want to save by setting `n_save`. Then, `n_save` files will be recorded during the propagation (`savename_1.mat`, `savename_2.mat`, ...). This enables to analyze the evolution of the pulse. `n_save` must be less or equal than `n_steps`. If `savename` is defined but not `n_save`, then only the result at the end of the propagation will be saved in a .mat file. This result is the parameter returned by the MEX function `mm_propagation`.

An example of how to set input parameters can be found in the **matlab** folder.
You can edit this file and adjust parameters for your simulation, or take it as a template.

# CONTACTS

MM_Propagation was developped by Alexandre Roget. You can contact him at alexandre.roget@inria.fr.
