# Modeling alpha-beta phase transformations in Ti-64

### Table of Contents  
1. [Introduction]()
2. [Scope]()
3. [Model - modules]()
4. [Usage instructions]()

### Introduction


### Scope


### Model - modules


### Usage instructions


### Using FFTW and This Code with Parallel

Things that need to be done:

* Need to convert the variants to the FFTW data structure.  There is an incompatibility between FFTW and MMSP.
* Two data structures are not compatible.
* Simplest working prototype using the SAME data structure:
  * Take a square wave
  * Diffusion step
  * Compute the FFT
  * Perform a filtering operations in the FFT space
  * Reverse the FFT
  * Diffusion step
* FFTW versus explicit calculation
