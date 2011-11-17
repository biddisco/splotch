# include "Galaxy.h"

float box_muller(float m, float s);

void CalculateColours (long npart, float * cred, float * cgreen, float * cblue, float * ciii,
                       float * Red, float * Green, float * Blue, float * III, float * xcoord, 
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

	      
	      cred[particlei]   = xcol;
	      cgreen[particlei] = xcol;
	      cblue[particlei]  = xcol;
	      ciii[particlei]   = xcol;

	   } else {
	      iaux = ii + jj*nx;
	      cred[particlei]   = Red[iaux];
	      cgreen[particlei] = Green[iaux];
	      cblue[particlei]  = Blue[iaux];
	      ciii[particlei]   = III[iaux];
	   }
	

	}
}
