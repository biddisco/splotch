#ifdef USEMPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <vector.h>


using namespace std;

/*
long Bin_reader (vector<float> &xpos, vector<float> &ypos, vector<float> &zpos,
                 vector<float> &scalar, vector<float> &smooth, float *maxr, float *minr,
                 long pe_size, int mype)
*/

long Bin_reader (vector<float> * xpos, vector<float> * ypos, vector<float> * zpos,
                 vector<float> * scalar, vector<float> * smooth, float *maxr, float *minr,
                 long pe_size, int mype)
{

   FILE * pFile;
   FILE * auxFile;
   float * dataarray;
   float * destarray;
   
   long total_size=0;
   float minradius=1e30;
   float maxradius=-1e30;

   int nfields = 6;
   float readarray[nfields];

   char datafile[500];
   long stride=pe_size*sizeof(float)*nfields*mype+4;

   scalar->resize(pe_size);
   xpos->resize(pe_size);
   ypos->resize(pe_size);
   zpos->resize(pe_size);
   smooth->resize(pe_size);


   if(mype == 0)
   { 
     printf("Input File Name\n");
     scanf("%s", datafile);
   }
#ifdef USEMPI
   MPI_Bcast(&datafile[0], 500, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

   cout << "DATAFILE INSIDE " << datafile << "     " << mype << "\n";

   pFile = fopen(datafile, "rb");

   cout << "Reading " << nfields << " fields for " << pe_size << " particles from " << stride << "\n";
   for(long index=0; index<pe_size; index++)
   {

      fseek(pFile, stride, SEEK_SET);
      fread(readarray, sizeof(readarray), 1, pFile);
      xpos->at(index) = readarray[0];
      ypos->at(index) = readarray[1];
      zpos->at(index) = readarray[2];
      scalar->at(index) = readarray[3];
      smooth->at(index) = readarray[5];
      //smooth.at(index) = 0.01;
      stride += sizeof(readarray);
      minradius = (minradius <= smooth->at(index) ? minradius : smooth->at(index));
      maxradius = (maxradius >= smooth->at(index) ? maxradius : smooth->at(index));

 
   }
   fclose(pFile);

   //maxradius = 1.0;
   *maxr=maxradius;
   *minr=minradius;

#ifdef USEMPI
   MPI_Allreduce(&maxradius, maxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&minradius, minr, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
   return pe_size;
}
