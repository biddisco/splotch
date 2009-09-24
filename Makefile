# type 'make ARCH=xxx'
# this is the basic makefile to compile splotch
# before using it, you have to create your Makefile.ARCH file in the makefiles directory
# where ARCH is a indentifier of your system. In Makefile.ARCH you have to define:
# SPLOTCH_DIR = splotch installation directory
# CXX = the compiler/loader
# OPT = optimization flags
# HDF5 = path to HDF5 library
# ZLIB = path to ZLIB library
# MPILIBS =  path to MPI libraries

include ./makefiles/Makefile.$(ARCH)

SPLOTCH_DIR = /sfs/sanfs/home/usercin/agh0/ParallelSplotch/splotch
INCLUDE = -I. -I$(SPLOTCH_DIR)/cxxsupport -I$(SPLOTCH_DIR)/kernel -I$(HDF5)/include
LIBS = -L$(HDF5)/lib -lhdf5 -L$(ZLIB)/lib -lz 
DEFINE = -D USEMPI -D SPLOTCH -DVAR_BRIGHT -D BCX -D DEBUG

EXEC = splotch$(ARCH)

$(EXEC) :
	$(CXX) $(OPT) $(DEFINE) $(INCLUDE) -o $(EXEC) fullsplotch.cc $(MPILIBS) $(LIBS)
