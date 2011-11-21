# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string>
# include <iostream>
# include "cxxsupport/paramfile.h"
# include "hdf5.h"

using namespace std;

long GaussRFunc (paramfile &params, string ComponentName, long number_of_points, long tot, float * coordx,
                 float * coordy, float * coordz,
                 float xmax, float ymax, float zmax, float * ddd, float * II, long nnx, long nny);

void CalculateDensity (float * hsml, float * rho, float * xcoord, float * ycoord,
                       float * zcoord, long numofpart, float smooth);

void CalculateColours (long npart, float * cred, float * cgreen, float * cblue, float * ciii, 
                       float * Red, float * Green, float * Blue, float * III, float * xcoord, 
                       float * ycoord, long nxxx, long nyyy);

long GlobularCluster (paramfile &params, string ComponentName, long number_of_points, long ntot,
                      float * coordx, float * coordy, float * coordz);

long ReadImages (paramfile &params, long numx, long numy, float Rmin, float * RRR,
                 float * GGG, float * BBB, float * III, float * xx, float * yy);

const int NUM_OF_FIELDS = 11;
