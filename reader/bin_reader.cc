#ifdef USE_MPI
#include "mpi.h"
#endif
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>

#ifdef VS
#include "cxxsupport/arr.h"
#include "cxxsupport/cxxutils.h"
#include "cxxsupport/paramfile.h"
#endif

#include "kernel/bstream.h"

using namespace std;
using namespace RAYPP;

#include "splotch/splotchutils.h"

long bin_reader_tab (paramfile &params, vector<particle_sim> &points, 
                     float *maxr, float *minr, 
                     int mype, int npes)
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
   long last_pe_adding;
   float minradius=1e30;
   float maxradius=-1e30;
// offset could be a input parameter
   long offset = 0;

   int n_load_fields = 8;
   int * which_fields = new int [n_load_fields];
   long totalsize;
   long totalsize_f;
   int num_of_fields;

   bool doswap = params.find<bool>("swap_endian",true);   
   string datafile = params.find<string>("infile");
   
   if(mype == 0)
   {
     num_of_fields = params.find<int>("num_columns");

     which_fields[0] = params.find<int>("x",-1);
     which_fields[1] = params.find<int>("y",-1);
     which_fields[2] = params.find<int>("z",-1);
     which_fields[3] = params.find<int>("r",-1);
     which_fields[4] = params.find<int>("I",-1);
     which_fields[5] = params.find<int>("C1",-1);
     which_fields[6] = params.find<int>("C2",-1);
     which_fields[7] = params.find<int>("C3",-1); 

     cout << "TABULAR BINARY FILE\n";
     cout << "Input data file name: " << datafile << endl;
     cout << "Number of columns " << num_of_fields << endl;
     cout << "x column (1 - " << num_of_fields << "), " << which_fields[0] << endl;
     cout << "y column (2 - " << num_of_fields << "), " << which_fields[1] << endl;
     cout << "z column (3 - " << num_of_fields << "), " << which_fields[2] << endl;
     cout << "r column (4 - " << num_of_fields << "), " << which_fields[3] << endl;
     cout << "I column (5 - " << num_of_fields << "), " << which_fields[4] << endl;
     cout << "C1 column (6 - " << num_of_fields << "), " << which_fields[5] << endl;
     cout << "C2 column (7 - " << num_of_fields << "), " << which_fields[6] << endl;
     cout << "C3 column (8 - " << num_of_fields << "), " << which_fields[7] << endl;

     for (int ii=0; ii<8; ii++)which_fields[ii]--;

     pFile = fopen(datafile.c_str(), "rb");
     fseek (pFile, 0, SEEK_END);
     totalsize = ftell (pFile);
     fclose(pFile);
   }
#ifdef USE_MPI
   MPI_Bcast(&num_of_fields, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(which_fields, n_load_fields, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&totalsize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
   // number of elements (of each variable) for a processor
   totalsize_f = (totalsize-offset)/(sizeof(float)*num_of_fields);
   pe_size = totalsize_f/npes;
   last_pe_adding = totalsize_f-pe_size*npes;
   if(mype == npes-1)pe_size += last_pe_adding;
   points.resize(pe_size);
   long pe_size_orig = pe_size;
#ifdef DEBUG
   if(mype == 0)  cout << "-----------------> " << pe_size << " " << last_pe_adding << "\n";
#endif
   long long stride=pe_size_orig*sizeof(float)*num_of_fields*mype+offset;

   float * readarray;
   readarray = new float [num_of_fields];
   bifstream infile;
   infile.open(datafile.c_str(), doswap);
   infile.rewind();
   
#ifdef DEBUG
   cout << "Reading " << num_of_fields << " fields for " << pe_size << " particles\n";
#endif

   infile.skip(stride);
   cout << mype << " -----> " << stride << "\n";
   for(long index=0; index<pe_size; index++)
   {

      for(int ifield=0; ifield<num_of_fields; ifield++) infile >> readarray[ifield];

      points.at(index).x=readarray[which_fields[0]];
      points.at(index).y=readarray[which_fields[1]];
      points.at(index).z=readarray[which_fields[2]];
      if(which_fields[3] >= 0)
      {
        points.at(index).r=readarray[which_fields[3]];
      } else {
        points.at(index).r=1.0;
      }
      if(which_fields[4] >= 0)
      {
        points.at(index).I=readarray[which_fields[4]];
      } else {
        points.at(index).I=0.5;
      }

      points.at(index).C1=readarray[which_fields[5]];

      if(which_fields[6] >= 0 && which_fields[7] >= 0)
      {
        points.at(index).C2=readarray[which_fields[6]];
        points.at(index).C3=readarray[which_fields[7]];
      } else {
        points.at(index).C2 = 0.0;
        points.at(index).C3 = 0.0;
      }
      points.at(index).ro=points.at(index).r;
      points.at(index).type=0;

      float smooth = points.at(index).r;      

      //stride += sizeof(float)*num_of_fields;
      minradius = (minradius <= smooth ? minradius : smooth);
      maxradius = (maxradius >= smooth ? maxradius : smooth);
 
   }
   infile.close();

   //if(mype == npes-1){for (int ii=0; ii<100000; ii+=1000)cout << mype << " " << points.at(ii).x << " " << points.at(ii).y << " " << points.at(ii).C1 << "\n";};
   //maxradius = 1.0;
   *maxr=maxradius;
   *minr=minradius;

#ifdef USE_MPI
   MPI_Allreduce(&maxradius, maxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&minradius, minr, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
   delete [] which_fields;
   delete [] readarray;
   return pe_size;
}

long bin_reader_block (paramfile &params, vector<particle_sim> &points, 
                       float *maxr, float *minr, 
                       int mype, int npes)
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
   int * which_fields = new int [n_load_fields];
   long totalsize;
   long totalsize_f;
   long last_pe_adding;

   bool doswap = params.find<bool>("swap_endian",true);   
   string datafile = params.find<string>("infile");
   int num_of_fields = params.find<int>("num_blocks",1);

   which_fields[0] = params.find<int>("x",-1);
   which_fields[1] = params.find<int>("y",-1);
   which_fields[2] = params.find<int>("z",-1);
   which_fields[3] = params.find<int>("r",-1);
   which_fields[4] = params.find<int>("I",-1);
   which_fields[5] = params.find<int>("C1",-1);
   which_fields[6] = params.find<int>("C2",-1);
   which_fields[7] = params.find<int>("C3",-1);

   if(mype == 0)
   {
     cout << "BLOCK BINARY FILE\n";
     cout << "Input data file name: " << datafile << endl;
     cout << "Number of blocks " << num_of_fields << endl;
     cout << "x block (1 - " << num_of_fields << "), " << which_fields[0] << endl;
     cout << "y block (2 - " << num_of_fields << "), " << which_fields[1] << endl;
     cout << "z block (3 - " << num_of_fields << "), " << which_fields[2] << endl;
     cout << "r block (4 - " << num_of_fields << "), " << which_fields[3] << endl;
     cout << "I block (5 - " << num_of_fields << "), " << which_fields[4] << endl;
     cout << "C1 block (6 - " << num_of_fields << "), " << which_fields[5] << endl;
     cout << "C2 block (7 - " << num_of_fields << "), " << which_fields[6] << endl;
     cout << "C3 block (8 - " << num_of_fields << "), " << which_fields[7] << endl;
     pFile = fopen(datafile.c_str(), "rb");
     fseek (pFile, 0, SEEK_END);
     totalsize = ftell (pFile);
     fclose(pFile);
     field_size = totalsize/(sizeof(float)*num_of_fields);
   }
#ifdef USE_MPI
   MPI_Bcast(&field_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

// number of elements (of each variable) for a processor
   pe_size = (long)(field_size / npes);
   last_pe_adding = field_size-pe_size*npes;
   long pe_size_orig = pe_size;
   if(mype == npes-1)pe_size += last_pe_adding;
   points.resize(pe_size);
   float * readarray;
   readarray = new float [pe_size];

   cout << "DATAFILE INSIDE " << datafile << "     " << mype << "\n";

   
   bifstream infile;
   infile.open(datafile.c_str(), doswap);

#ifdef DEBUG
   cout << "Reading " << n_load_fields << " fields for " << pe_size << " particles\n";
#endif
   for(int n_fields=0; n_fields<n_load_fields; n_fields++)
   {
     int n_fields_eff = which_fields[n_fields]-1;
     if(which_fields[n_fields] < 0)continue;

     stride=sizeof(float)*(n_fields_eff*field_size+pe_size_orig*mype)+offset;

     infile.rewind();
     infile.skip(stride);
     for(long index=0; index<pe_size; index++)infile >> readarray[index];

     switch(n_fields)
     {
     case 0:
       for(long index=0; index<pe_size; index++)points.at(index).x=readarray[index];
       break;
     case 1:
       for(long index=0; index<pe_size; index++)points.at(index).y=readarray[index];
       break;
     case 2:
       for(long index=0; index<pe_size; index++)points.at(index).z=readarray[index];
       break;
     case 3:
       for(long index=0; index<pe_size; index++)points.at(index).r=readarray[index];
       break;
     case 4:
       for(long index=0; index<pe_size; index++)points.at(index).I=readarray[index];
       break;
     case 5:
       for(long index=0; index<pe_size; index++)points.at(index).C1=readarray[index];
       break;
     case 6:
       for(long index=0; index<pe_size; index++)points.at(index).C2=readarray[index];
       break;
     case 7:
       for(long index=0; index<pe_size; index++)points.at(index).C3=readarray[index];
       break;
     }
     for(long index=0; index<pe_size; index++)
     {
       points.at(index).ro=points.at(index).r;
       points.at(index).type=0;
       float smooth = points.at(index).r;      

       minradius = (minradius <= smooth ? minradius : smooth);
       maxradius = (maxradius >= smooth ? maxradius : smooth);
     }

   }
   infile.close();


   if(which_fields[4] < 0)
       for(long index=0; index<pe_size; index++)points.at(index).I=0.5;

   //maxradius = 1.0;
   *maxr=maxradius;
   *minr=minradius;

#ifdef USE_MPI
   MPI_Allreduce(&maxradius, maxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&minradius, minr, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif
   delete [] which_fields;
   delete [] readarray;
   return pe_size;
}


