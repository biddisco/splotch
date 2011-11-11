# include <algorithm>
# include "Galaxy.h"


long ReadImages (paramfile &params, long numx, long numy, unsigned int Rmin, unsigned int * RRR,
                 unsigned int * GGG, unsigned int * BBB, float * xx, float * yy)
{

	  FILE * pFile;
	  float toll = 0.0;
          unsigned char * RR;
          unsigned char * GG;
          unsigned char * BB;
          float * x;
	  float * y;
	  unsigned int RRaux;

// Read Data: this part will be changed with FITSIO reading three different files
// one for each color. For testing, we just refer to the single raw file

          string infile = params.find<string>("GalaxyFileR","NONE");
	  pFile = fopen(infile.c_str(), "rb");

// image size will be read from the FITS file: at the moment it's an input

//       read numx
//       read numy
	 long num = numx*numy;
	 RR = new unsigned char[num];
	 GG = new unsigned char[num];
	 BB = new unsigned char[num];
	 x = new float [num];
	 y = new float [num];

	 long counter = 0;
	 float dist = 0.0;
	 float xaux, yaux;
	 float xxx=0.0;
	 float yyy=0.0;

// set the center of the galaxy (at the moment in term of pixels (default center of image)
         float norm0 = (float)numx;
	 float xc0 = (float)numx/2/norm0;
	 float yc0 = (float)numy/2/norm0;
         float norm = params.find<float>("Gnorm", norm0);
         float xc = params.find<float>("Gcenterx",xc0);
         float yc = params.find<float>("Gcentery",yc0);

         long i=0;
          for (int iy=0; iy<numy; iy++)
          for (int ix=0; ix<numx; ix++)
	  { 
		fread(&RR[i], sizeof(char), 1, pFile);
		fread(&GG[i], sizeof(char), 1, pFile);
		fread(&BB[i], sizeof(char), 1, pFile);
                i++;
          }                 
	  fclose(pFile);

// When FITSIO implemented other two files to be read

          infile = params.find<string>("GalaxyFileG","NONE");
// Read Green file

          infile = params.find<string>("GalaxyFileB","NONE");
// Read Blue file

// Now data processing

          i=0;
          int THR = 0;
          for (int iy=0; iy<numy; iy++)
          for (int ix=0; ix<numx; ix++)
          {
		RRR[i] = (unsigned int)RR[i];
		GGG[i] = (unsigned int)GG[i];
		BBB[i] = (unsigned int)BB[i];

		x[i] = (float)ix;
		y[i] = (float)iy;
		x[i] = x[i]/norm-xc;
		y[i] = y[i]/norm-yc;

                THR = max((int)RRR[i], (int)GGG[i]);
                THR = max(THR, (int)BBB[i]);
		if(THR >= Rmin && RRR[i] < 257)
///////		if(RRR[i] >= Rmin)
		{
		   xx[counter] = x[i];
		   yy[counter] = y[i];
		   counter++;
		}
		i++;
	  }

	  delete [] RR;
	  delete [] GG;
	  delete [] BB;

	  printf("=== NUMBER OF ACTIVE PIXELS : %ld ===\n", counter);
	  return (counter);

}
