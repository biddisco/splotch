#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <vector.h>

#include "arr.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"

using namespace std;
using namespace RAYPP;

#include "splotch/splotchutils.h"

long bin_reader_tab (vector<particle_sim> &points, float *maxr, float *minr, 
                     char*datafile, int num_of_fields, 
                     int * which_fields, int mype, int npes)
{
/*
In this case we expect the file to be written as a binary table
xyzrabcde
xyzrabcde
xyzrabcde
...
xyzrabcde

which_fields represents the column for the various quantity according
to the following standard:
which_fields[0] = x coord
which_fields[1] = y coord
which_fields[2] = z coord
which_fields[3] = r coord
which_fields[4] = intensity
which_fields[5] = color 1 (R)
which_fields[6] = color 2 (G)
which_fields[7] = color 3 (B)
*/


   FILE * pFile;
   FILE * auxFile;
   float * dataarray;
   float * destarray;
   
   long total_size=0;
   long pe_size;
   float minradius=1e30;
   float maxradius=-1e30;
// offset could be a input parameter
   long offset = 0;

   int n_load_fields = 8;
   long totalsize;
   if(mype == 0)
   {
     pFile = fopen(datafile, "rb");
     fseek (pFile, 0, SEEK_END);
     totalsize = ftell (pFile);
     fclose(pFile);
// number of elements (of each variable) for a processor
     pe_size = totalsize/(sizeof(float)*npes*num_of_fields);
   }
#ifdef USE_MPI
   MPI_Bcast(&pe_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

   points->resize(pe_size);
   float * readarray;
   readarray = new float [num_of_fields];

   long stride=pe_size*sizeof(float)*num_of_fields*mype+offset;

#ifdef USE_MPI
   MPI_Bcast(&datafile[0], 500, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

   cout << "DATAFILE INSIDE " << datafile << "     " << mype << "\n";

   pFile = fopen(datafile, "rb");

#ifdef DEBUG
   cout << "Reading " << nfields << " fields for " << pe_size << " particles from " << stride << "\n";
#endif
   for(long index=0; index<pe_size; index++)
   {

      fseek(pFile, stride, SEEK_SET);
      fread(readarray, sizeof(readarray), 1, pFile);
      points->at(index).x=readarray[which_fields[0]];
      points->at(index).y=readarray[which_fields[1]];
      points->at(index).z=readarray[which_fields[2]];
      points->at(index).r=readarray[which_fields[3]];
      points->at(index).I=readarray[which_fields[4]];
      points->at(index).C1=readarray[which_fields[5]];
      points->at(index).C2=readarray[which_fields[6]];
      points->at(index).C3=readarray[which_fields[7]];
      points->at(index).ro=points->at(index).r;

      float smooth = points->at(index).r;      

      stride += sizeof(readarray);
      minradius = (minradius <= smooth ? minradius : smooth);
      maxradius = (maxradius >= smooth ? maxradius : smooth);
 
   }
   fclose(pFile);

   //maxradius = 1.0;
   *maxr=maxradius;
   *minr=minradius;

#ifdef USE_MPI
   MPI_Allreduce(&maxradius, maxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&minradius, minr, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
   delete [] readarray;
   return pe_size;
}

long bin_reader_block (vector<particle_sim> &points, float *maxr, float *minr, 
                     char*datafile, int num_of_fields, 
                     int * which_fields, int mype, int npes)
{
/*
In this case we expect the file to be written as 
xxxxxxxxxx 
xxxxxxxxxx 
yyyyyyyyyy
yyyyyyyyyy
...
TTTTTTTTTT
TTTTTTTTTT

which_fields represents the block position inside the file according
to the following standard:
which_fields[0] = x coord
which_fields[1] = y coord
which_fields[2] = z coord
which_fields[3] = r coord
which_fields[4] = intensity
which_fields[5] = color 1 (R)
which_fields[6] = color 2 (G)
which_fields[7] = color 3 (B)
*/

   FILE * pFile;
   FILE * auxFile;
   float * dataarray;
   float * destarray;
   
   long total_size=0;
   long pe_size;
   long field_size;
   float minradius=1e30;
   float maxradius=-1e30;
// offset could be a input parameter
   long offset = 0;
   long stride;

   int n_load_fields = 8;
   long totalsize;
   if(mype == 0)
   {
     pFile = fopen(datafile, "rb");
     fseek (pFile, 0, SEEK_END);
     totalsize = ftell (pFile);
     fclose(pFile);
// number of elements (of each variable) for a processor
     field_size = totalsize/(sizeof(float)*num_of_fields);
     pe_size = (long)(field_size / npes);
   }
#ifdef USE_MPI
   MPI_Bcast(&pe_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

   points->resize(pe_size);
   float * readarray;
   readarray = new float [pe_size];

#ifdef USE_MPI
   MPI_Bcast(&datafile[0], 500, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

   cout << "DATAFILE INSIDE " << datafile << "     " << mype << "\n";

   pFile = fopen(datafile, "rb");

#ifdef DEBUG
   cout << "Reading " << nfields << " fields for " << pe_size << " particles from " << stride << "\n";
#endif
   for(int n_fields=0; n_fields<n_load_fields; n_fields++)
   {
     int n_fields_eff = which_fields[n_fields];
     stride=sizeof(float)*(n_fields_eff*field_size+pe_size*mype)+offset;

     fseek(pFile, stride, SEEK_SET);
     fread(readarray, sizeof(readarray), 1, pFile);
     switch(n_fields)
     {
     case 0:
       for(long index=0; index<pe_size; index++)points->at(index).x=readarray[index];
       break;
     case 1:
       for(long index=0; index<pe_size; index++)points->at(index).y=readarray[index];
       break;
     case 2:
       for(long index=0; index<pe_size; index++)points->at(index).z=readarray[index];
       break;
     case 3:
       for(long index=0; index<pe_size; index++)points->at(index).r=readarray[index];
       break;
     case 4:
       for(long index=0; index<pe_size; index++)points->at(index).I=readarray[index];
       break;
     case 5:
       for(long index=0; index<pe_size; index++)points->at(index).C1=readarray[index];
       break;
     case 6:
       for(long index=0; index<pe_size; index++)points->at(index).C2=readarray[index];
       break;
     case 7:
       for(long index=0; index<pe_size; index++)points->at(index).C3=readarray[index];
       break;
     }
     for(long index=0; index<pe_size; index++)
     {
       points->at(index).ro=points->at(index).r;
       float smooth = points->at(index).r;      

       minradius = (minradius <= smooth ? minradius : smooth);
       maxradius = (maxradius >= smooth ? maxradius : smooth);
     }

   }
   fclose(pFile);

   //maxradius = 1.0;
   *maxr=maxradius;
   *minr=minradius;

#ifdef USE_MPI
   MPI_Allreduce(&maxradius, maxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&minradius, minr, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
   delete [] readarray;
   return pe_size;
}


