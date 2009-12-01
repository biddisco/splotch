/*
 ============================================================================
 Name        : partition.c
 Author      : GFM
 Version     : 0.01
 Copyright   :
 Description : MPI2 I/O
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "partition.h"


int
MPI_Ajo_partition (int nprocs, int nsize[], int npxd[], int npart)
{

    int i, j;
    int *x;

    x = (int *) malloc (nprocs * sizeof (int));

    npxd[2] = 1;
// partizionamento di una sola dimensione
    if (npart == 1)
    {
        npxd[0] = nprocs;
        for (i = 1; i < 2; ++i)
            npxd[i] = 1;
    }
// partizionamento di due dimensioni
    else if (npart == 2)
    {

        for (i = 1; i <= nprocs; ++i)
        {
            x[i - 1] = 0;
            if (nprocs % i == 0)
                x[i - 1] = i;
        }

        int piccolo = nprocs + 2;

        for (i = 0; i < nprocs; ++i)
        {
            for (j = i; j < nprocs; ++j)
            {
                if (x[i] * x[j] == nprocs)
                {
                    if (x[i] + x[j] <= piccolo)
                    {
                        npxd[0] = x[i];
                        npxd[1] = x[j];
                        piccolo = x[i] + x[j];
                    }
                }
            }
        }

        int temp;
        if (nsize[0] >= nsize[1])
        {
            temp = npxd[0];
            npxd[0] = npxd[1];
            npxd[1] = temp;
        }
    }
    else
    {
        printf ("Inserire un numero corretto (1 o 2)\n");
        return 1;
    }

    free (x);
    return 0;
}


int
MPI_Ajo_new_array (MPI_Comm comm, int ndim_array, int nsize_global[],
                   int nsize[], int psize[], int my_rank,
                   int start_global_array[])
{
    int i, j;
    int *nresto;
    int *indice;
    int ierr;
    int nsize_local = 0;
    int *periods;
    int reorder = 0;

    MPI_Comm comm_cart;

    periods = (int *) malloc (ndim_array * sizeof (int));
    indice = (int *) malloc (ndim_array * sizeof (int));
    nresto = (int *) malloc (ndim_array * sizeof (int));

    for (i = 0; i < ndim_array; ++i)
        periods[i] = 0;

    ierr =
        MPI_Cart_create (comm, ndim_array, psize, periods, reorder, &comm_cart);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Cart_create");

    ierr = MPI_Cart_coords (comm_cart, my_rank, ndim_array, indice);
    MPI_Ajo_msgerr (comm_cart, my_rank, ierr, "MPI_Cart_coords");

//      printf ( "INDICI=%d %d %d %d\n",indice[0],indice[1],indice[2],my_rank );

    for (i = 0; i < ndim_array; ++i)
    {
        start_global_array[i] = 0;
        nresto[i] = nsize_global[i] % psize[i];

        nsize[i] = nsize_global[i] / psize[i];
        if (indice[i] < nresto[i])
            ++nsize[i];

        for (j = 0; j < indice[i]; ++j)
        {
            nsize_local = nsize_global[i] / psize[i];
            if (j < nresto[i])
                ++nsize_local;
            start_global_array[i] += nsize_local;
        }
    }

    free (periods);
    free (indice);
    free (nresto);
    return 0;
}

int
MPI_Ajo_write (MPI_Comm comm, int my_rank, const char *filename, int ndim_array,
               int nsize_global[], int nsize[], int start_global_array[],
               MPI_Datatype etype, void *array, MPI_Offset offset)
{
    int ierr = 0;

    MPI_Datatype filetype;
    MPI_Status status;
    MPI_File fh;

    int i, count = 1;
    for (i = 0; i < ndim_array; ++i)
        count *= nsize[i];


    ierr =
        MPI_File_open (comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                       MPI_INFO_NULL, &fh);
       MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_open");
    ierr =
        MPI_Type_create_subarray (ndim_array, nsize_global, nsize,
                                  start_global_array, MPI_ORDER_C, etype,
                                  &filetype);
#ifdef DEBUG
       MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Type_create_subarray");

    ierr =
        MPI_Ajo_check_size (comm, my_rank, &etype, &filetype, ndim_array, nsize);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Ajo_check_size");
#endif
    ierr = MPI_Type_commit (&filetype);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Type_commit");

    ierr = MPI_File_seek (fh, offset, MPI_SEEK_SET);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_seek");

    ierr =
        MPI_File_set_view (fh, offset, etype, filetype, "native", MPI_INFO_NULL);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_set_view");

    ierr = MPI_File_write (fh, array, count, etype, &status);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_write");

    ierr = MPI_File_close (&fh);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_close");
    return ierr;
}

int
MPI_Ajo_read (MPI_Comm comm, int my_rank, const char *filename, int ndim_array,
              int nsize_global[], int nsize[], int start_global_array[],
              MPI_Datatype etype, void *array, MPI_Offset offset)
{
    int ierr = 0;

    MPI_Datatype filetype;
    MPI_Status status;
    MPI_File fh;

    int i, count = 1;
    for (i = 0; i < ndim_array; ++i)
        count *= nsize[i];

    ierr = MPI_File_open (comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
   MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_open");

    ierr =
        MPI_Type_create_subarray (ndim_array, nsize_global, nsize,
                                  start_global_array, MPI_ORDER_C, etype,
                                  &filetype);
#ifdef DEBUG
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Type_create_subarray");

    ierr =
        MPI_Ajo_check_size (comm, my_rank, &etype, &filetype, ndim_array, nsize);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Ajo_check_size");
#endif
    ierr = MPI_Type_commit (&filetype);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Type_commit");

    ierr = MPI_File_seek (fh, offset, MPI_SEEK_SET);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_seek");

    ierr =
        MPI_File_set_view (fh, offset, etype, filetype, "native", MPI_INFO_NULL);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_set_view");

    ierr = MPI_File_read (fh, array, count, etype, &status);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_read");

    ierr = MPI_File_close (&fh);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_File_close");

    return MPI_SUCCESS;
}

int
MPI_Ajo_check_size (MPI_Comm comm, int my_rank, MPI_Datatype * etype,
                    MPI_Datatype * filetype, int ndim_array, int nsize[])
{
    int i;
    int ierr = 0;

    int size_etype;
    int size_filetype;
    int size_filetype_calc = 1;

    ierr = MPI_Type_size (*filetype, &size_filetype);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Type_size");

    ierr = MPI_Type_size (*etype, &size_etype);
    MPI_Ajo_msgerr (comm, my_rank, ierr, "MPI_Type_size");

    for (i = 0; i < ndim_array; ++i)
        size_filetype_calc *= nsize[i];
    size_filetype_calc *= size_etype;

    if (size_filetype_calc != size_filetype)
        return -1;

    return MPI_SUCCESS;
}

void
my_MPI_Ajo_msgerr (MPI_Comm comm, int my_rank, int ierr, char *stringa,
                   const char *filename, const char *funct_name, int nline)
{
#ifdef DEBUG
    int total_errors = 0;

    MPI_Allreduce (&ierr, &total_errors, 1, MPI_INT, MPI_SUM, comm);

    if (total_errors != MPI_SUCCESS)
    {

        if (my_rank == 0)
            fprintf (stderr, "ERROR in file %s, function %s, line %d, by %s\n",
                     filename, funct_name, nline, stringa);
        fflush (stderr);

        MPI_Finalize ();
        exit (0);
    }
    else
        return;
#endif
}
