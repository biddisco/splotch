#######################################################################
#  Splotch V4                                                         #
#######################################################################

#--------------------------------------- Basic operation mode of code
#OPT	+=  -DGEOMETRY_FILE
#OPT	+=  -DINTERPOLATE
#OPT	+=  -DHIGH_ORDER_INTERPOLATION


#--------------------------------------- Switch on MPI
OPT	+=  -DUSE_MPI
OPT	+=  -DUSE_MPIIO

#--------------------------------------- Visual Studio Option
#OPT	+=  -DVS

#--------------------------------------- CUDA options
#OPT	+=  -DCUDA
#OPT	+=  -DCUDA_THREADS
#OPT     +=  -DNO_WIN_THREAD
#OPT     +=  -DNO_HOST_RANGING
#OPT     +=  -DNO_HOST_TRANSFORM
#OPT     +=  -DNO_HOST_COLORING
#OPT     +=  -DNO_HOST_RENDER
#OPT	 +=  -DHOST_THREAD_RENDER
#OPT	 +=  -DCUDA_DEVICE_COMBINE
#OPT	 +=  -DCUDA_TEST_COLORMAP
#OPT	 +=  -DCUDA_TEST_FRAGMENT


#--------------------------------------- Select target Computer

SYSTYPE="SP6"
#SYSTYPE="GP"
#SYSTYPE="PLX"

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       = mpic++        # sets the C-compiler (default)
else
CC       = g++        # sets the C-compiler (default)
endif
OMP      = -fopenmp

OPTIMIZE = -Wextra -Wall -Wstrict-aliasing=2 -Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wold-style-cast -O2 -g    # optimization and warning flags (default)
SUP_INCL = -I. -Icxxsupport -Impiio-1.0/include/

ifeq ($(SYSTYPE),"SP6")
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  mpCC_r
else
CC       =  xlc++        
endif
OPTIMIZE =  -q64 -O3 -qarch=auto -qtune=auto -qinline
LIB_OPT	 =  -bstackpsize:64k -bdatapsize:64k -btextpsize:64k
OMP =
endif

ifeq ($(SYSTYPE),"GP")
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  /opt/cuda/bin/nvcc
else
CC       =  /opt/cuda/bin/nvcc
endif
OPTIMIZE = -O2 
LIB_OPT  = -Xlinker -L -Xlinker /opt/cuda/lib
OMP =  
#-Xcompiler -openmp
SUP_INCL += -I/opt/cuda/sdk/common/inc -I/opt/cuda/include -Icuda
endif

ifeq ($(SYSTYPE),"PLX")
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  nvcc -g
else
CC       =  nvcc -g
endif
OPTIMIZE = -O2 -DDEBUG
LIB_OPT  = -Xlinker -L$(NVCC_HOME)/lib
OMP =
SUP_INCL += -I$(CUDASDK_HOME)/common/inc -I$(NVCC_HOME)/include -Icuda
endif

#--------------------------------------- Here we go

OPTIONS = $(OPTIMIZE) $(OPT)

EXEC   = Splotch4.0$(SYSTYPE)

OBJS  =	kernel/transform.o utils/colourmap.o cxxsupport/error_handling.o \
	cxxsupport/mpi_support.o cxxsupport/cxxutils.o reader/gadget_reader.o \
	reader/millenium_reader.o reader/bin_reader.o reader/bin_reader_mpi.o \
	writer/write_tga.o splotch/splotchutils.o splotch/splotch.o

ifeq (CUDA,$(findstring CUDA,$(OPT)))
OBJS += cuda/splotch.o cuda/CuPolicy.o
endif
ifeq (USE_MPIIO,$(findstring USE_MPIIO,$(OPT)))
LIB_MPIIO = -Lmpiio-1.0/lib -lpartition
endif

INCL   = splotch/splotchutils.h writer/writer.h reader/reader.h	Makefile

CPPFLAGS = $(OPTIONS) $(SUP_INCL) $(OMP)

CUFLAGS = $(OPTIONS) $(SUP_INCL)

LIBS   = $(LIB_OPT) $(OMP)

.SUFFIXES: .o .cc .cxx .cpp .cu

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cpp.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cu.o:
	$(CC) -c $(CUFLAGS) -o "$@" "$<"

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC) $(LIB_MPIIO)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)

