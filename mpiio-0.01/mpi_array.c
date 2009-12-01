/*
 ============================================================================
 Name        : mpi_array.c
 Author      : GFM
 Version     : 0.01
 Copyright   :
 Description : MPI2 I/O
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "partition.h"

int main ( int argc, char* argv[] )
{

	int i,j,k;
	int ndim_array;
	int *nsize_global;
	int tot_dati_global;
	int *start_global_array;
	int *psize;

	int *nsize;
	float *array_monodim;
	int tot_dati;
	float ***array;

	int ierr=0;
	int my_rank=0;          /* rank of process */
	int nprocs;           /* number of processes */

	MPI_Status status ;   
	MPI_Datatype etype=MPI_FLOAT;	
	MPI_Offset   offset=0;

	char *mpi_file_write, *mpi_file_write2;
	mpi_file_write="write1.bin";
	mpi_file_write2="write2.bin";

/*START MPI*/
	ierr = MPI_Init ( &argc, &argv );

	ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Comm_rank" );

	ierr = MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Comm_size" );

  
 /*LETTURA FILE INPUT - START*/

	MPI_File fparam;

	ierr = MPI_File_open ( MPI_COMM_WORLD, "bin.param",MPI_MODE_RDONLY, MPI_INFO_NULL, &fparam );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_open" );

	ierr = MPI_File_seek ( fparam, 0*sizeof ( int ), MPI_SEEK_SET );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_seek" );

	ierr = MPI_File_read ( fparam,&ndim_array,1,MPI_INT,&status );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_read" );

	nsize              = ( int * ) malloc ( ndim_array * sizeof ( int ) );
	nsize_global       = ( int * ) malloc ( ndim_array * sizeof ( int ) );
	psize              = ( int * ) malloc ( ndim_array * sizeof ( int ) );
	start_global_array = ( int * ) malloc ( ndim_array * sizeof ( int ) );

	for ( i=0;i<ndim_array; ++i )
	{
		ierr = MPI_File_seek ( fparam, ( 1+i ) *sizeof ( int ), MPI_SEEK_SET );
		MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_seek" );

		ierr = MPI_File_read ( fparam,&nsize_global[i],1,MPI_INT,&status );
		MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_read" );
	}
	ierr = MPI_File_close ( &fparam );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_close" );
/*LETTURA FILE INPUT - END*/
// // // // // // // // // //


	printf ( "<<GLOBAL: %d %d %d %d %d>>\n", ndim_array, nsize_global[0], nsize_global[1], nsize_global[2], my_rank );

/*FUNZIONE DI PARTIZIONAMENTO*/
	ierr = MPI_Ajo_partition ( nprocs, nsize_global, psize, 2 );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Ajo_partition" );

	ierr = MPI_Ajo_new_array ( MPI_COMM_WORLD, ndim_array, nsize_global, nsize, psize, my_rank , start_global_array );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Ajo_new_array" );

	printf ( "<<LOCAL  : %d %d %d %d %d>>\n", ndim_array,        nsize[0],        nsize[1],        nsize[2], my_rank );
	printf ( "<<GLOBAL : %d %d %d %d %d>>\n", ndim_array, start_global_array[0], start_global_array[1], start_global_array[2], my_rank );

/*FUNZIONE DI PARTIZIONAMENTO - END*/

/*ALLOC DELLA MATRICE LOCALE*/
	tot_dati = nsize[0]*nsize[1]*nsize[2];
	tot_dati_global = nsize_global[0]*nsize_global[1]*nsize_global[2];
	array_monodim = ( float * ) malloc ( tot_dati * sizeof ( float ) );
	array = ( float *** ) malloc ( nsize[0] * sizeof ( float ** ) );
	for ( i=0;i<nsize[0];++i )
		array[i] = ( float ** ) malloc ( nsize[1] * sizeof ( float * ) );
	for ( i=0;i<nsize[0];++i )
		for ( j=0;j<nsize[1];++j )
			array[i][j] = &array_monodim[nsize[1]*nsize[2]*i+j*nsize[2]];

/*INIZIALIZZAZIONE DELLA MATRICE LOCALE*/
	for ( i=0;i<nsize[0];++i )
		for ( j=0;j<nsize[1];++j )
			for ( k=0;k<nsize[2];++k )
				array[i][j][k]=my_rank;



	ierr = MPI_Ajo_write ( MPI_COMM_WORLD, my_rank, mpi_file_write, ndim_array, nsize_global, nsize, start_global_array, etype, &array[0][0][0], offset );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Ajo_write" );

	offset=tot_dati_global*sizeof ( float );

	ierr = MPI_Ajo_write ( MPI_COMM_WORLD, my_rank, mpi_file_write, ndim_array, nsize_global, nsize, start_global_array, etype, &array[0][0][0], offset );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Ajo_write" );


/*AZZERAMENTO MATRICE*/
	offset=0;
	for ( i=0;i<nsize[0];++i )
		for ( j=0;j<nsize[1];++j )
			for ( k=0;k<nsize[2];++k )
			{
				array[i][j][k]=0;
			}

/*LETTURA FILE*/
	ierr = MPI_Ajo_read ( MPI_COMM_WORLD, my_rank, mpi_file_write, ndim_array, nsize_global, nsize, start_global_array, etype, &array[0][0][0], offset );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_File_read" );

/*SCRITTURA SU ALTRO FILE*/
	ierr = MPI_Ajo_write ( MPI_COMM_WORLD, my_rank, mpi_file_write2, ndim_array, nsize_global, nsize, start_global_array, etype, &array[0][0][0], offset );
	MPI_Ajo_msgerr ( MPI_COMM_WORLD , my_rank, ierr, "MPI_Ajo_write" );

	for ( i=0;i<nsize[0];++i )
		free ( array[i] );
	free ( array );

	free ( array_monodim );
	free ( nsize );
	free ( nsize_global );
	free ( psize );
	free ( start_global_array );

	MPI_Finalize();

	return 0;
}
