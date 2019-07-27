# Makefile

# includes
MMSP_PATH=/home/arun/mmsp
incdir = $(MMSP_PATH)/include
utildir = $(MMSP_PATH)/utility
algodir = $(MMSP_PATH)/algorithms

# compilers/flags
compiler = g++ -std=c++11 -O3 
pcompiler = mpic++ -O3
flags = -I$(incdir) -I$(algodir) -I$(utildir)
prefourierflags = -I/home/arun/Documents/mmsp-arun/fftw-3.3.8/api -fopenmp 
postfourierflags = /home/arun/Documents/mmsp-arun/fftw++-2.05/fftw++.cc  
fourierlinkers = -lfftw3 -lfftw3_omp -lm

# IBM compiler for AMOS
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_XL)/include
BG_LIB = -L$(BG_XL)/lib
qcompiler = $(BG_XL)/bin/mpixlcxx_r -O3 -qflag=w -qstrict -qmaxmem=-1
qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib

# dependencies
core = #$Thermo.hpp \
       #$(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.h \
       $(incdir)/MMSP.grid.h \
       $(incdir)/MMSP.sparse.h \

# the program
sparse: ti64_master.cpp $(core)
	$(pcompiler) $(flags) $(fourierflags) -include mpi.h $(fourierlinkers) $< -o sparse.out -lz

sparse_serial: $(core)
	$(compiler) $(flags) $(prefourierflags) ti64_anisose.cpp $(postfourierflags) $(fourierlinkers) $< -o sparse_serial.out -lz

bgq: sp-xmpf.cpp $(core)
	$(qcompiler) -DBGQ $(qflags) $< -o q_sparse.out -lz

clean:
	rm -rf sparse_serial.out pf_* delta* c_* *~
