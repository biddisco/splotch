#######################################################################
#  Splotch V5                                                      #
#######################################################################

#--------------------------------------- Switch on DataSize
#OPT += -DLONGIDS

#--------------------------------------- Switch on MPI
OPT += -DUSE_MPI
#OPT += -DUSE_MPIIO

#--------------------------------------- Switch on HDF5

#OPT += -DHDF5
#OPT += -DH5_USE_16_API

#--------------------------------------- Visual Studio Option
#OPT += -DVS

#--------------------------------------- CUDA options

OPT += -DCUDA
OPT += -DNO_WIN_THREAD

#--------------------------------------- OpenCL options

#OPT += -DOPENCL
#OPT += -DNO_WIN_THREAD

#--------------------------------------- Turn on VisIVO stuff
#OPT += -DSPLVISIVO

#--------------------------------------- Turn off Intensity  normalization
#OPT += -DNO_I_NORM

#--------------------------------------- Select target Computer
#SYSTYPE="generic"
#SYSTYPE="SP6"
#SYSTYPE="GP"
SYSTYPE="PLX"
#SYSTYPE="BGP"
#SYSTYPE="VIZ"
#SYSTYPE="EIGER"
#SYSTYPE="TODI"
### visualization cluster at the Garching computing center (RZG):
#SYSTYPE="RZG-SLES11-VIZ"
### generic SLES11 Linux machines at the Garching computing center (RZG):
#SYSTYPE="RZG-SLES11-generic"




# Set compiler executables to commonly used names, may be altered below!
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
 CC       = mpic++
else
 CC       = g++
endif

# optimization and warning flags (g++)
OPTIMIZE = -std=c++98 -pedantic -Wno-long-long -Wfatal-errors -Wextra -Wall -Wstrict-aliasing=2 -Wundef -Wshadow -Wwrite-strings -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wold-style-cast
ifeq ($(SYSTYPE),"generic")
# OPTIMIZE += -O2 -g -D TWOGALAXIES
 OPTIMIZE += -O2 -g
endif

# OpenMP compiler switch
OMP      = #-fopenmp

SUP_INCL = -I. -Icxxsupport -Ic_utils

#CUDA_HOME = /usr/local/cuda/
ifeq (USE_MPIIO,$(findstring USE_MPIIO,$(OPT)))
 SUP_INCL += -Impiio-1.0/include/
endif


# Configuration for the VIZ visualization cluster at the Garching computing centre (RZG):
# ->  gcc/OpenMPI_1.4.2
ifeq ($(SYSTYPE),"RZG-SLES11-VIZ")
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC       = mpic++
 else
  CC       = g++
 endif
 ifeq (HDF5,$(findstring HDF5,$(OPT)))
  HDF5_HOME = /u/system/hdf5/1.8.7/serial
  LIB_HDF5  = -L$(HDF5_HOME)/lib -Wl,-rpath,$(HDF5_HOME)/lib -lhdf5 -lz
  HDF5_INCL = -I$(HDF5_HOME)/include
 endif
 OPTIMIZE += -O3 -march=native -mtune=native
 OMP       = -fopenmp
endif

# configuration for TODI at CSCS
ifeq ($(SYSTYPE),"TODI")
 ifeq (HDF5,$(findstring HDF5,$(OPT)))
  HDF5_HOME = /opt/cray/hdf5/1.8.6/gnu/46/
  LIB_HDF5  = -L$(HDF5_HOME)/lib -lhdf5 -lz
  HDF5_INCL = -I$(HDF5_HOME)/include
 endif
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
   CC       = CC 
 else
   CC       = g++
 endif
 OPTIMIZE = -O3
 OMP      = 
endif

# Configuration for SLES11 Linux clusters at the Garching computing centre (RZG):
# ->  gcc/IntelMPI_4.0.0, requires "module load impi"
ifeq ($(SYSTYPE),"RZG-SLES11-generic")
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC       = mpigxx
 else
  CC       = g++
 endif
 ifeq (HDF5,$(findstring HDF5,$(OPT)))
  HDF5_HOME = /afs/ipp/home/k/khr/soft/amd64_sles11/opt/hdf5/1.8.7
  LIB_HDF5  = -L$(HDF5_HOME)/lib -Wl,-rpath,$(HDF5_HOME)/lib -lhdf5 -lz
  HDF5_INCL = -I$(HDF5_HOME)/include
 endif
 OPTIMIZE += -O3 -msse3
 OMP       = -fopenmp
endif


ifeq ($(SYSTYPE),"SP6")
 ifeq (HDF5,$(findstring HDF5,$(OPT)))
  HDF5_HOME = /cineca/prodDF5_INCL = -I$(HDF5_HOME)/include/libraries/hdf5/1.8.4_ser/xl--10.1
  LIB_HDF5  = -L$(HDF5_HOME)/lib -lhdf5 -L/cineca/prod/libraries/zlib/1.2.3/xl--10.1/lib/ -lz -L/cineca/prod/libraries/szlib/2.1/xl--10.1/lib/ -lsz
 endif
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC       =  mpCC_r
 else
  CC       =  xlc++
 endif
 OPTIMIZE =  -q64 -O3 -qarch=auto -qtune=auto -qinline
 LIB_OPT	 =  -bstackpsize:64k -bdatapsize:64k -btextpsize:64k
 OMP =
endif

ifeq ($(SYSTYPE),"EIGER")
ifeq (HDF5,$(findstring HDF5,$(OPT)))
HDF5_HOME = /scratch/eiger/cgheller/hdf5-1.8.7/install/
LIB_HDF5  = -L$(HDF5_HOME)/lib -lhdf5
HDF5_INCL = -I$(HDF5_HOME)/include
endif
ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
CC       =  mpic++
else
CC       =  c++
endif
#EIGER again
ifeq (CUDA,$(findstring CUDA,$(OPT)))
NVCC       =  nvcc -g
CUDA_HOME  =  /apps/eiger/Cuda-4.0/cuda
#LIB_OPT  += -L$(CUDA_HOME)/lib -L$(CUDA_HOME)/lib64 -lcudart
LIB_OPT  += -L$(CUDA_HOME)/lib64 -lcudart
SUP_INCL += -I/apps/eiger/NVIDIA_GPU_Computing_SDK/C/common/inc/ -I$(CUDA_HOME)/include#-I$(CUDAUTIL_INC) -I$(NVCC_HOME)/include
endif

OPTIMIZE =  -O3
OMP =
endif

ifeq ($(SYSTYPE),"GP")
 CC       =  nvcc -g
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC       =  mpicxx -g -I$(CUDA_HOME)/sdk/common/inc -I$(CUDA_HOME)/sdk/C/common/inc -I$(CUDA_HOME)/include
 endif
 NVCC       =  nvcc -g 
 OPTIMIZE = -O2
 LIB_OPT  =  -L$(CUDA_HOME)/lib -L$(CUDA_HOME)/lib64 -lcudart
 OMP =  
 #-Xcompiler -openmp
 SUP_INCL += -I$(CUDA_HOME)/sdk/common/inc -I$(CUDA_HOME)/sdk/C/common/inc # -I$(CUDA_HOME)/include  -Icuda
endif

ifeq ($(SYSTYPE),"BGP")
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC       = mpixlcxx_r
 else
  CC       = bgxlC_r
 endif
 OPTIMIZE = -O3 -qstrict -qarch=450d -qtune=450 # -qipa=inline
 LIB_OPT  =
 OMP =   -qsmp=omp -qthreaded
endif

ifeq ($(SYSTYPE),"PLX")
 ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g 
 else
  CC  = g++
 endif
 OPTIMIZE = -O2 -DDEBUG
 OMP = #-fopenmp
 ifeq (CUDA,$(findstring CUDA,$(OPT)))
  NVCC = nvcc -arch sm_20 -use_fast_math
  LIB_OPT  =  -L$(CUDA_HOME)/lib64 -lcudart
  SUP_INCL += -I$(CUDA_HOME)/include -I$(CUDA_SDK)/CUDALibraries/common/inc 
 endif
 ifeq (OPENCL,$(findstring OPENCL,$(OPT)))
  LIB_OPT  =  -L$(CUDA_HOME)/lib64 -L$(CUDA_HOME)/lib -lOpenCL
  SUP_INCL += -I$(CUDA_HOME)/include
 endif
endif

#-L/home/pavel/NVIDIA_GPU_Computing_SDK/shared/lib 
#
#--------------------------------------- Here we go

OPTIONS = $(OPTIMIZE) $(OPT)

EXEC = Splotch5-$(SYSTYPE)
EXEC1 = Galaxy

OBJS  =	kernel/transform.o cxxsupport/error_handling.o \
        reader/mesh_reader.o reader/visivo_reader.o \
	cxxsupport/mpi_support.o cxxsupport/paramfile.o cxxsupport/string_utils.o \
	cxxsupport/announce.o cxxsupport/ls_image.o reader/gadget_reader.o \
	reader/millenium_reader.o reader/bin_reader.o reader/bin_reader_mpi.o reader/tipsy_reader.o \
	splotch/splotchutils.o splotch/splotch.o \
	splotch/scenemaker.o splotch/splotch_host.o cxxsupport/walltimer.o c_utils/walltime_c.o \
	booster/mesh_creator.o booster/randomizer.o booster/p_selector.o booster/m_rotation.o

OBJS1 = galaxy/Galaxy.o galaxy/GaussRFunc.o galaxy/Box_Muller.o galaxy/ReadBMP.o \
	galaxy/CalculateDensity.o galaxy/CalculateColours.o galaxy/GlobularCluster.o \
	galaxy/ReadImages.o

OBJSC = cxxsupport/paramfile.o cxxsupport/error_handling.o cxxsupport/mpi_support.o \
	c_utils/walltime_c.o cxxsupport/string_utils.o \
	cxxsupport/announce.o cxxsupport/ls_image.o cxxsupport/walltimer.o

ifeq (HDF5,$(findstring HDF5,$(OPT)))
 OBJS += reader/hdf5_reader.o
 OBJS += reader/gadget_hdf5_reader.o
 OBJS += reader/galaxy_reader.o
endif

ifeq (OPENCL,$(findstring OPENCL,$(OPT)))
 OBJS += opencl/splotch.o opencl/CuPolicy.o opencl/splotch_cuda2.o opencl/deviceQuery.o
else
 ifeq (CUDA,$(findstring CUDA,$(OPT)))
  OBJS += cuda/splotch.o cuda/CuPolicy.o cuda/splotch_cuda2.o cuda/deviceQuery.o
 endif
endif

ifeq (USE_MPIIO,$(findstring USE_MPIIO,$(OPT)))
 LIB_MPIIO = -Lmpiio-1.0/lib -lpartition
endif

INCL   = */*.h Makefile

CPPFLAGS = $(OPTIONS) $(SUP_INCL) $(HDF5_INCL) $(OMP)

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
	$(NVCC) -c $(CUFLAGS) -o "$@" "$<"

$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC) $(LIB_MPIIO) $(LIB_HDF5)

$(EXEC1): $(OBJS1) $(OBJSC)
	$(CC) $(OPTIONS) $(OBJS1) $(OBJSC) $(LIBS) -o $(EXEC1) $(LIB_HDF5)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS)

cleangalaxy:
	rm -f $(OBJS1)

realclean: clean
	rm -f $(EXEC)

