#include "kernel/transform.cxx"
#include "utils/colourmap.cxx"
#include "cxxsupport/error_handling.cc"
#include "cxxsupport/mpi_support.cc"
#include "cxxsupport/cxxutils.cc"
#include "reader/gadget_reader.cc"
#include "reader/millenium_reader.cc"
//#include "reader/enzo_reader.cc"
#include "reader/bin_reader.cc"
#ifdef USE_MPI
#include "reader/bin_reader_mpi.cc"
#endif
#include "writer/write_tga.cc"
#include "splotch/splotch.cc"
