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
* Parallel: Runs without fftw. The roadblock here is to sync the way mmsp handles parallelism and fftw handles parallelism. 


### Model - modules


### Usage instructions

Please refer to the compilation guide. 
