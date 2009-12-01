AC_DEFUN([AC_ENABLE_DISABLE],
[
#####################################################################################
#			  ENABLE    /   DISABLE
#####################################################################################
 
#####################################################################################
## ENABLE/DISABLE DEBUG
AC_ARG_ENABLE(debug,AC_HELP_STRING([--enable-debug=yes|no],[Enable support for debugging (default=no)]),
  [enable_debug=$enableval],
  [enable_debug=no]
)
if test "$enable_debug" = "yes"; then
AC_MSG_NOTICE([===============])
AC_MSG_NOTICE([  ENABLE DEBUG])
AC_MSG_NOTICE([===============])
CFLAGS=" -O0 -g -Wall"
CXXFLAGS=" -O0 -g -Wall"
FFLAGS=" -O0 -g"
FCFLAGS=" -O0 -g"
fi
AM_CONDITIONAL([DEBUG], [test $enable_debug = yes])

#####################################################################################
## ENABLE TYPE
AC_ARG_ENABLE(type,AC_HELP_STRING([--enable-type=serial|seek_wr|wr_at|setv_wr|darray_wr],[Enable choise of version (default=serial)]),
  [enable_type=$enableval],
  [enable_type="serial darray_wr"]
)
serial=0
mpi=0
seek_wr=0
wr_at=0
setv_wr=0
darray_wr=0

for gg in $enable_type;do
	case "$gg" in
	serial)
		serial=yes ;;
	seek_wr)
		seek_wr=yes ;;
	wr_at)
		wr_at=yes ;;
	setv_wr)
		setv_wr=yes ;;
	darray_wr)
		darray_wr=yes ;;
	esac
done


if test "$serial" = "yes"; then
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([ENABLE SERIAL VERSION])
AC_MSG_NOTICE([=====================])
fi
AM_CONDITIONAL([SERIAL], [test $serial = yes])

if test "$seek_wr" = "yes"; then
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([ENABLE SEEK + WRITE VERSION])
AC_MSG_NOTICE([=====================])
AC_MPI
CC=$MPICC
mpi=yes
fi

if test "$wr_at" = "yes"; then
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([ENABLE WRITE_AT VERSION])
AC_MSG_NOTICE([=====================])
AC_MPI
CC=$MPICC
mpi=yes
fi

if test "$setv_wr" = "yes"; then
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([ENABLE SETVIEW + WRITE VERSION])
AC_MSG_NOTICE([=====================])
AC_MPI
CC=$MPICC
mpi=yes
fi

if test "$darray_wr" = "yes"; then
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([ENABLE DARRAY + WRITE VERSION])
AC_MSG_NOTICE([=====================])
AC_MPI
CC=$MPICC
mpi=yes
fi

if test "$mpi" = "yes"; then
AC_MSG_NOTICE([=====================])
AC_MSG_NOTICE([ENABLE MPI VERSION])
AC_MSG_NOTICE([=====================])
fi

AM_CONDITIONAL([SEEK_WR], [test $seek_wr = yes ])
AM_CONDITIONAL([WR_AT], [test $wr_at = yes ])
AM_CONDITIONAL([SETV_WR], [test $setv_wr = yes ])
AM_CONDITIONAL([DARRAY_WR], [test $darray_wr = yes ])
AM_CONDITIONAL([MPI], [test $mpi = yes])


]) 