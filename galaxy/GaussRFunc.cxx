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

  printf("      disk with sigma_z = %f, %d particles per pixel\n", sigma_z, n_per_pixel);

  long count = 0;

  float* coordx_save;
  float* coordy_save;
  coordx_save = new float [number_of_points];
  coordy_save = new float [number_of_points];

  for (long i=0; i<number_of_points; i++)
    {
      coordx_save[i] = coordx[i];
      coordy_save[i] = coordy[i];
    }

  float pixsizex = 2./nx;
  float pixsizey = 2./ny;
  float coordz_aux;
  for (long i=0; i<number_of_points; i++)
    {
      float x0=coordx_save[i];
      float y0=coordy_save[i];

      coordx[count] = x0;
      coordy[count] = y0;
      coordz[count] = box_muller(0.0, sigma_z);
      coordz_aux = coordz[count];
      
      count++;
      if(count >= ntot)
	{
	  printf("Generating more particles than allowed (%d>=%d)\n",count,ntot);
	  exit(3);
	}

      for(long k=1; k<n_per_pixel; k++)
	{
	  //coordx[count] = box_uniform(x0, pixsizex);
	  //coordy[count] = box_uniform(y0, pixsizey);
	  coordx[count] = box_muller(x0, pixsizex);
	  coordy[count] = box_muller(y0, pixsizey);
	  // the factor of 4.0 just below is set "experimetally"
	  coordz[count] = box_muller(coordz_aux, sigma_z/4.0);
	  count++;
	  if(count >= ntot)
	    {
	      printf("Generating more particles than allowed (%d>=%d)\n",count,ntot);
	      exit(3);
	    }
	}
    }

  delete [] coordx_save;
  delete [] coordy_save;

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
	  if(III[rindex] > 0.0)
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

  float xmin=1e10;
  float xmax=-1e10;
  float ymin=1e10;
  float ymax=-1e10;
  for (long i=0; i<number_of_points; i++)    // copy pixel positions to the leading of the array
    {
      if(coordz[i] == 0.0)
	continue;
      coordx[iaux] = coordx[i];
      coordy[iaux] = coordy[i];
      coordz[iaux] = coordz[i];
      xmax = max(coordx[iaux],xmax);
      ymax = max(coordy[iaux],ymax);
      xmin = min(coordx[iaux],xmin);
      ymin = min(coordy[iaux],ymin);
      iaux++;
    }

  float ref_size = ((xmax-xmin)+(ymax-ymin))/2.0;
  cout << "CHARCHTERISTIC SIZE = " << ref_size << "\n";

  long abscounter = countin;

  long ii = 0;
  float pixsizex = 4./rx;
  float pixsizey = 4./ry;
  while (abscounter < ntot-npergroup && ii < countin)
    {
      //original by Klaus
      //float gsigmaz = sigma * (sigma_fixed + coordy[ii]);
      //Claudio:
      float weight = coordz[ii];
      float thick = sqrt(weight)*ref_size;
      //float gsigmaz = sqrt(sigma) * 0.25* thick;
      float gsigmaz = sigma * thick;
      

      coordz[ii] = box_muller(0.0, gsigmaz);
      for(long jj=1; jj<(long)(weight*npergroup); jj++)
	{
	  //coordx[abscounter] = box_uniform(coordx[ii], pixsizex);
	  //coordy[abscounter] = box_uniform(coordy[ii], pixsizey);
	  coordx[abscounter] = box_muller(coordx[ii], pixsizex*1.0);
	  coordy[abscounter] = box_muller(coordy[ii], pixsizey*1.0);
	  float zaux;
	  do {zaux = box_muller(coordz[ii], sigma_fixed*gsigmaz);}
	  //////////do {zaux = box_muller(coordz[ii], sigma_fixed);}
	  //while (zaux*zaux > 0.25*ref_size*ref_size);
	  while (zaux*zaux > 0.65*ref_size*ref_size);
          coordz[abscounter] = zaux;

	  abscounter++;
	}
      ii++;
    }

  return abscounter;

}
