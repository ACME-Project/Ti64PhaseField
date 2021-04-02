# Makefile

# includes
MMSP_PATH= ./mmsp
incdir = $(MMSP_PATH)/include
utildir = $(MMSP_PATH)/utility
algodir = $(MMSP_PATH)/algorithms

# compilers/flags
compiler = g++ -std=c++11 -O3 
pcompiler = mpic++ -std=c++11 -O3
flags = -I$(incdir) -I$(algodir) -I$(utildir)
prefourierflags = -I./fftw-3.3.9/api -fopenmp 
postfourierflags = ./fftw++-2.05/fftw++.cc  
fourierlinkers = -lfftw3 -lfftw3_omp -lm  
fourier_mpi_linkers = -lfftw3 -lfftw3_mpi -lm  

# dependencies
core = #$Thermo.hpp \
       #$(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.h \
       $(incdir)/MMSP.grid.h \
       $(incdir)/MMSP.sparse.h \

# Serial compilation - without fftw. Remove the fftw-related functions in the code before you can use this compilation 
#sparse: $(core)
#	$(compiler) ./Serial/main.cpp $< -o sparse.out -lz  #-include mpi.h

#Serial compilation - with fftw
sparse_serial: $(core)
	$(compiler) $(flags) $(prefourierflags) ./Serial/main.cpp $(postfourierflags) $(fourierlinkers)  $< -o sparse_serial.out -lz #
	
#Parallel compilation - without fftw
parallel: $(core)
	$(pcompiler) $(flags) ./Parallel/main.cpp $< -o parallel.out -lz  #-include mpi.h
	
test: $(core)
	$(pcompiler) $(flags) $(prefourierflags)  ./Parallel_fftw/main.cpp $(postfourierflags) $(fourier_mpi_linkers) $< -o test.out -lz  #-include mpi.h



bgq: sp-xmpf.cpp $(core)
	$(qcompiler) -DBGQ $(qflags) $< -o q_sparse.out -lz

tool: tool.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(pcompiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT

mmsp2png : mmsp2png.cpp
	$(compiler) $(flags) $< -o $@ -lz -lpng

# convert MMSP grid file to ParaView Data file type
mmsp2pvd : mmsp2pvd.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to tab-delimited ASCII (TSV) file type
mmsp2tsv : mmsp2tsv.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to VTK Image file type
mmsp2vti : mmsp2vti.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to XYZ point cloud file type
mmsp2xyz : mmsp2xyz.cpp
	$(compiler) $(flags) $< -o $@ -lz

clean:
	rm -rf sparse.out output_* pf_* delta* c_* IE* *~ variant*
	
	

#Parallel compilation - with fftw (not recommended. Integration of mmsp with fftw on a parallel environment has not been tested)
#parallel: $(core)
	#$(pcompiler) $(flags) $(prefourierflags)  "Add your path for main.cpp here" $(postfourierflags) $(fourierlinkers) $< -o parallel.out -lz  #-include mpi.h
