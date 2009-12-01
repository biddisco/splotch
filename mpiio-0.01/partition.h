/*
 ============================================================================
 Name        : partition.c
 Author      : GFM
 Version     : 0.01
 Copyright   :
 Description : MPI2 I/O
 ============================================================================
 */

#ifndef partition_
#define partition_

#ifdef __cplusplus
extern "C"
{
#endif

#define MPI_Ajo_msgerr(x,y,z,zz) my_MPI_Ajo_msgerr(x,y,z,zz, __FILE__, __FUNCTION__, __LINE__);


int MPI_Ajo_partition ( int nprocs, int nsize[], int npxd[], int npart );

int MPI_Ajo_new_array ( MPI_Comm comm, int ndim_array, int nsize_global[], int nsize[], int psize[], int my_rank, int start_global_array[] );

int MPI_Ajo_write ( MPI_Comm comm, int my_rank, const char *filename, int ndim_array, int nsize_global[], int nsize[], int start_global_array[], MPI_Datatype etype, void *array, MPI_Offset offset );

int MPI_Ajo_read ( MPI_Comm comm, int my_rank, const char *filename, int ndim_array, int nsize_global[], int nsize[], int start_global_array[], MPI_Datatype etype, void *array, MPI_Offset offset );

int MPI_Ajo_check_size ( MPI_Comm comm, int my_rank, MPI_Datatype *etype, MPI_Datatype *filetype, int ndim_array, int nsize[] );

void my_MPI_Ajo_msgerr ( MPI_Comm comm, int my_rank, int ierr, char *stringa ,const char *filename, const char *funct_name, int nline );


#ifdef __cplusplus
}
#endif

#endif /*partition_*/

