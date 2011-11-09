# include "Galaxy.h"

float box_muller(float m, float s);

void CalculateColours (long npart, unsigned int * cred, unsigned int * cgreen, unsigned int * cblue, unsigned int * ciii,
                       unsigned int * Red, unsigned int * Green, unsigned int * Blue, unsigned int * III, float * xcoord, 
                       float * ycoord, long nx, long ny)
{

	float xaux, yaux, xcol;
	long ii, jj, iaux;
        float x_rand_max = (float) RAND_MAX;
	float xcolaux;

	for (long particlei=0; particlei<npart; particlei++)
	{

	   xaux = (0.5*(xcoord[particlei]+1.0)); 
	   yaux = (0.5*(ycoord[particlei]+1.0)); 
	   ii = (int) (xaux*nx);
	   jj = (int) (yaux*ny);

	   if(ii >= nx || ii < 0 || jj >=ny || jj < 0)
	   {
//              xcol = ((float)rand())/x_rand_max;
	      xcolaux = box_muller(0, 0.25);
	      xcol = fabs(xcolaux);
	      if (xcol > 1.0) xcol = 0.0;

	      
	      cred[particlei]   = (unsigned int)(xcol*255);
	      cgreen[particlei] = (unsigned int)(xcol*255);
	      cblue[particlei]  = (unsigned int)(xcol*255);
	      ciii[particlei]   = (unsigned int)(xcol*255);

	   } else {
	      iaux = ii + jj*nx;
	      cred[particlei]   = Red[iaux];
	      cgreen[particlei] = Green[iaux];
	      cblue[particlei]  = Blue[iaux];
	      ciii[particlei]   = III[iaux];
	   }
	

	}
}\
