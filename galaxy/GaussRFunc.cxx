# include "Galaxy.h"

float box_muller(float m, float s);
float box_uniform(float m, float s);

long GaussRFunc (paramfile &params, string ComponentName, long number_of_points, 
		 float * coordx, float * coordy, float * coordz) 
{
  float sigma[3];
	
  srand(time(NULL));

  sigma[0] = params.find<float>("Sigmax"+ComponentName,0);
  sigma[1] = params.find<float>("Sigmay"+ComponentName,0);
  sigma[2] = params.find<float>("Sigmaz"+ComponentName,0);

  printf("    spheroid component with sigma_[x,y,z] = %f,%f,%f\n", sigma[0],sigma[1],sigma[2]);

// x/y/z coord

  for (long i=0; i<number_of_points; i++)
    {
      coordx[i] = box_muller(0.0, sigma[0]);
      coordy[i] = box_muller(0.0, sigma[1]);
      coordz[i] = box_muller(0.0, sigma[2]);
    }

  return number_of_points;

}

long GaussRDiscFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz, long nx, long ny) 
{
  srand(time(NULL));

  float sigma_z = params.find<float>("Sigmaz"+ComponentName,0.01);
  long n_per_pixel = params.find<float>("NperPixel"+ComponentName,1);

  printf("      disk with sigma_z = %f\n", sigma_z);

  long count = 0;

// x/y/z coord

  float pixsizex = 2./nx;
  float pixsizey = 2./ny;
  for (long i=0; i<number_of_points; i++)
    {
      for(long k=0; k<n_per_pixel; k++)
	{
	  coordx[count] = box_uniform(coordx[i], pixsizex);
	  coordy[count] = box_uniform(coordy[i], pixsizey);
	  coordz[count] = box_muller(0.0, sigma_z);
	  count++;
	  if(count >= ntot)
	    {
	      printf("Generating more particles than allowed (%d>=%d)\n",count,ntot);
	      exit(3);
	    }
	}
    }

  return count;

}



long RDiscFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz,
                 float * III, long nx, long ny) 
{
  srand(time(NULL));

  float sigma = params.find<float>("Sigmaz"+ComponentName,0);
  float sigma_fixed = params.find<float>("Sigmazfixed"+ComponentName,0.1);
  long npergroup = params.find<long>("NperGroup"+ComponentName,0);
  long rx = params.find<long>("Scaledxres"+ComponentName,1);
  long ry = params.find<long>("Scaledyres"+ComponentName,1);
  float compression = params.find<float>("CompressionFactor"+ComponentName,1.0);
  float * resolution;

  resolution = new float [rx*ry];
  long countin = 0;
  long countout = 0;

  for(long i=0; i<rx*ry; i++)
    resolution[i] = 0;

  for (long i=0; i<number_of_points; i++)    // find unique points related to low res image description
    {
      int ix = (int) round((0.5*(coordx[i]+1.0))*(rx));
      int iy = (int) round((0.5*(coordy[i]+1.0))*(ry));
      int irx = (int) round((0.5*(coordx[i]+1.0))*nx);
      int iry = (int) round((0.5*(coordy[i]+1.0))*ny);

      long index = ix+iy*rx;
      long rindex = irx+iry*nx;
      
      if(resolution[index] == 0.0)
	{
	  if(III[rindex] > 0)
	    {
	      coordz[i] = III[rindex];
	      resolution[index] = 1.0;
	      countin++;
	    }
	  else
	    {
	      coordz[i] = 0.0;
	      countout++;
	    }
	} 
      else 
	{
	  coordz[i] = 0.0;
	  countout++;
	}
    }

  delete [] resolution;

  long iaux = 0;

  for (long i=0; i<number_of_points; i++)    // copy pixel positions to the leading of the array
    {
      if(coordz[i] == 0.0)
	continue;
      coordx[iaux] = coordx[i];
      coordy[iaux] = coordy[i];
      coordz[iaux] = coordz[i];
      iaux++;
    }

  long abscounter = countin;

  long ii = 0;
  float pixsizex = 2./rx;
  float pixsizey = 2./ry;
  while (abscounter < ntot-npergroup && ii < countin)
    {
      float gsigmaz = sigma * (sigma_fixed + coordy[ii]);
      coordz[ii] = box_muller(0.0, gsigmaz);
      for(long jj=0; jj<npergroup; jj++)
	{
	  coordx[abscounter] = box_uniform(coordx[ii], pixsizex);
	  coordy[abscounter] = box_uniform(coordy[ii], pixsizey);
	  coordz[abscounter] = box_muller(coordz[ii], gsigmaz);
	  abscounter++;
	}
      ii++;
    }

  return abscounter-1;

}
