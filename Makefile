# Makefile

MMSP_PATH = /home/arun/mmsp

# includes
incdir = $(MMSP_PATH)/include
utildir = $(MMSP_PATH)/utility
algodir = $(MMSP_PATH)/algorithms

# compilers/flags
compiler = g++ -O3 -Wall
pcompiler = mpic++ -O3
flags = -I$(incdir) -I$(algodir) -I$(utildir)

# IBM compiler for AMOS
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_XL)/include
BG_LIB = -L$(BG_XL)/lib
qcompiler = $(BG_XL)/bin/mpixlcxx_r -O3 -qflag=w -qstrict -qmaxmem=-1
qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.h \
       $(incdir)/MMSP.grid.h \
       $(incdir)/MMSP.sparse.h \

# the program
sparse: sp-xmpf.cpp sp-initialize.hpp $(core)
	$(pcompiler) $(flags) -include mpi.h $< -o sparse.out -lz

bgq: sp-xmpf.cpp $(core)
	$(qcompiler) -DBGQ $(qflags) $< -o q_sparse.out -lz

clean:
	rm -rf q_GG.out sparse.out *~
