#######################################################################
#  Splotch V4                                                         #
#######################################################################

#--------------------------------------- Basic operation mode of code
OPT	+=  -DGEOMETRY_FILE
OPT	+=  -DINTERPOLATE
OPT	+=  -DHIGH_ORDER_INTERPOLATION


#--------------------------------------- Switch on MPI
#OPT	+=  -DUSE_MPI
#OPT	+=  -DUSE_MPIIO

#--------------------------------------- Visual Studio Option
#OPT	+=  -DVS

#--------------------------------------- CUDA options
#OPT	+=  -DCUDA
#OPT	+=  -DHOST_THREAD_RENDER
#OPT	+=  -DCUDA_DEVICE_COMBINE
#OPT	+=  -DCUDA_THREADS
#OPT	+=  -DCUDA_TEST_COLORMAP
#OPT	+=  -DCUDA_TEST_FRAGMENT


#--------------------------------------- Select target Computer

#SYSTYPE="SP6"

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       = mpic++        # sets the C-compiler (default)
else
CC       = g++        # sets the C-compiler (default)
endif
OMP      = -fopenmp

OPTIMIZE = -Wall  -g    # optimization and warning flags (default)
SUP_INCL = -I. -Icxxsupport

ifeq ($(SYSTYPE),"SP6")
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  mpCC_r
else
CC       =  xlc++        
endif
OPTIMIZE =  -q64 -O3 -qarch=auto -qtune=auto -qinline
LIB_OPT	 =  -bstackpsize:64k -bdatapsize:64k -btextpsize:64k
endif

#--------------------------------------- Here we go

OPTIONS = $(OPTIMIZE) $(OPT)

EXEC   = Splotch4.0

OBJS  =	kernel/transform.o utils/colourmap.o cxxsupport/error_handling.o \
	cxxsupport/mpi_support.o cxxsupport/cxxutils.o reader/gadget_reader.o \
	reader/millenium_reader.o reader/bin_reader.o reader/bin_reader_mpi.o \
	writer/write_tga.o splotch/splotchutils.o splotch/splotch.o

INCL   = splotch/splotchutils.h writer/writer.h reader/reader.h	Makefile

CPPFLAGS = $(OPTIONS) $(SUP_INCL) $(OMP)

LIBS   = $(LIB_OPT) $(OPTIONS) $(OMP)

.SUFFIXES: .o .cc .cxx

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)

