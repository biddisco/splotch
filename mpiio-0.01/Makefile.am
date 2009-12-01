ACLOCAL_AMFLAGS = -I m4
m4dir = m4
include m4/Makefile.mk

bin_PROGRAMS = gen_param mpi_array
dist_bin_SCRIPTS = od.ksh
# od_ksh_SOURCES = od.ksh
# main_par main_partition main_ser mpi_array_1d

# if SERIAL
# endif

# if MPI
# endif

# main_ser_SOURCES = $(main_par_SOURCES)
# main_ser_CPPFLAGS = $(enable_debug)
# main_ser_LDADD =
# main_ser_CFLAGS =


# main_par_SOURCES = main_par.c
# main_par_CPPFLAGS = $(enable_mpi) $(enable_debug) $(enable_io)
# main_par_LDADD =
# main_par_CFLAGS =

mpi_array_SOURCES = mpi_array.c
mpi_array_CPPFLAGS = $(enable_mpi) $(enable_debug) $(enable_io)
mpi_array_LDADD = $(top_builddir)/libpartition.a
# CFLAGS = -g -O0 -Wall

lib_LIBRARIES = libpartition.a
libpartition_a_SOURCES = partition.c partition.h
include_HEADERS = partition.h

gen_param_SOURCES = asci.param gen_param.c
dist_doc_DATA     = asci.param
# docdir= $(dir)
# mpi_array_1d_SOURCES = MPI_array_1d.c
# mpi_array_1d_CPPFLAGS = $(enable_mpi) $(enable_debug) $(enable_io)

# main_partition_SOURCES = main_partition.c


# main_partition_LDADD = $(top_builddir)/libpartition.a
