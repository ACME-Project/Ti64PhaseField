# Modeling alpha-beta phase transformations in Ti-64

### Table of Contents  
1. [Config_params]()
2. [Tested compilations]()
3. [Model - modules]()
4. [Usage instructions]()

### Config_params

A config file, which will feed domain parameters, will be added in the next release. Until then, please use the definitions_v2.hpp file, lines 1-37, to modify domain parameters. Any discrepacy in the following parameters can especially lead to an error:

* dim
* nx, ny, and nz: Make sure that the value of dim is compatible with grid definition.

### Tested compilations

The following repos have been tested for compilation errors:

* Serial: fftw is included.
* Parallel: Runs without fftw. 


### Compilations under development

Parallel computation with fftw. 

The current stage of testing has completed the following steps:

1.) double and complex arrays have been defined under the rules of fftw parallelization
2.) dft plans have been defined for forward and inverse transforms under mpi environment

Next steps:

Reproduce the serial fftw computations.


### Model - modules


### Usage instructions

Please refer to the compilation guide. 
