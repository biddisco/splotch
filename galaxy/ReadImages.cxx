# include <algorithm>
# include "Galaxy.h"


long ReadImages (paramfile &params, long numx, long numy, float Rmin, float * RRR,
                 float * GGG, float * BBB, float * III, float * xx, float * yy)
{

	  FILE * pFile;
	  FILE * pFile1;
	  float toll = 0.0;
          unsigned char * RR;
          unsigned char * GG;
          unsigned char * BB;
          unsigned char * II;
          float * x;
	  float * y;

// Read Data: this part will be changed with FITSIO reading three different files
// one for each color. For testing, we just refer to the single raw file

          string infile = params.find<string>("GalaxyFileR","NONE");
	  pFile = fopen(infile.c_str(), "rb");
          string infile1 = params.find<string>("GalaxyFileI","NONE");
	  pFile1 = fopen(infile1.c_str(), "rb");

// image size will be read from the FITS file: at the moment it's an input

//       read numx
//       read numy
	 long num = numx*numy;
	 RR = new unsigned char[num];
	 GG = new unsigned char[num];
	 BB = new unsigned char[num];
         II = new unsigned char[num];
	 x = new float [num];
	 y = new float [num];

	 long counter = 0;
	 float dist = 0.0;
	 float xaux, yaux;
	 float xxx=0.0;
	 float yyy=0.0;

// set the center of the galaxy (at the moment in term of pixels (default center of image)
         float norm0 = 0.5*(float)numx;
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
		fread(&II[i], sizeof(char), 1, pFile1);
                i++;
          }                 
	  fclose(pFile);
	  fclose(pFile1);

// When FITSIO implemented other two files to be read

          infile = params.find<string>("GalaxyFileG","NONE");
// Read Green file

          infile = params.find<string>("GalaxyFileB","NONE");
// Read Blue file

// Now data processing

          i=0;
          float THR = 0;
          for (int iy=0; iy<numy; iy++)
          for (int ix=0; ix<numx; ix++)
          {

		RRR[i] = (float)RR[i]/255.0;
		GGG[i] = (float)GG[i]/255.0;
		BBB[i] = (float)BB[i]/255.0;
		III[i] = (float)II[i]/255.0;

		x[i] = (float)ix;
		y[i] = (float)iy;
		x[i] = x[i]/norm-xc;
		y[i] = y[i]/norm-yc;

                //THR = max(RRR[i], GGG[i]);
                //THR = max(THR, BBB[i]);
                THR = III[i];
 
		if(THR >= Rmin)
///////		if(THR >= Rmin && RRR[i] < 1.0)
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
	  delete [] II;

	  printf("=== NUMBER OF ACTIVE PIXELS : %ld ===\n", counter);

	  /*
          pFile = fopen("test.txt", "w");
          for(long ii=0; ii<counter; ii++)
          {
             fprintf(pFile, "%f %f 0.0\n", xx[ii],yy[ii]);
          }
          fclose(pFile);
	  */
	  return (counter);

}
