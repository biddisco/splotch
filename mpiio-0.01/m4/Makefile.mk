if SERIAL
enable_serial=-DSERIAL
endif

if MPI
enable_mpi=-DMPI
endif

if DEBUG
enable_debug=-DDEBUG
endif

if SEEK_WR
enable_io=-DSEEK_WR
endif

if WR_AT
enable_io=-DWR_AT
endif

if SETV_WR
enable_io=-DSETV_WR
endif

if DARRAY_WR
enable_io=-DDARRAY_WR
endif