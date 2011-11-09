# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string>
# include <iostream>
# include "cxxsupport/paramfile.h"
# include "hdf5.h"

using namespace std;


long ReadBMP (paramfile &params, string infile, string infile1, long numx, long numy, unsigned int Rmin,
              unsigned int Gmin, unsigned int Bmin, unsigned int * RR,
              unsigned int * GG, unsigned int * BB, unsigned int * II, float * xx, float * yy);

long GaussRFunc (paramfile &params, string ComponentName, long number_of_points, long tot, float * coordx,
                 float * coordy, float * coordz,
                 float xmax, float ymax, float zmax, float * ddd, unsigned int * II, long nnx, long nny);

void CalculateDensity (float * hsml, float * rho, float * xcoord, float * ycoord,
                       float * zcoord, long numofpart, float smooth);

void CalculateColours (long npart, unsigned int * cred, unsigned int * cgreen, unsigned int * cblue, unsigned int * ciii, 
                       unsigned int * Red, unsigned int * Green, unsigned int * Blue, unsigned int * III, float * xcoord, 
                       float * ycoord, long nxxx, long nyyy);

long GlobularCluster (paramfile &params, string ComponentName, long number_of_points, long ntot,
                      float * coordx, float * coordy, float * coordz);


const int NUM_OF_FIELDS = 11;

