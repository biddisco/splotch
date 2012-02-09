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

long GaussRGlobFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz, float * III, long nx, long ny) 
{
  srand(time(NULL));
  float * xx;
  float * yy;
  float * zz;

  long npergroup = params.find<long>("NperGroup"+ComponentName,0);
  float sigma_fixed = params.find<float>("Sigmazfixed"+ComponentName,0.1);


// set the center of the galaxy (at the moment in term of pixels (default center of image)
  float norm = 0.5*(float)nx;
  float xc = (float)nx/2/norm;
  float yc = (float)ny/2/norm;

  FILE * pFile;
// read local maxima in the pixels distribution
/*
	string imagefile_max = params.find<string>(ComponentName+"FileMax","NONE");
        unsigned char * max_mask = new unsigned char[nx*ny];
        pFile = fopen(imagefile_max.c_str(), "rb");
        int maxcount = 0;
        for (int ix=0; ix<nx; ix++)
        for (int iy=0; iy<ny; iy++)
           {
             long index = ix+iy*nx;
             fread(&max_mask[index], sizeof(char), 1, pFile);
             //if (max_mask[index] != 0)maxcount++;
             if(max_mask[index] > 100)maxcount++;
           }
        fclose(pFile);
*/

// find local maxima in the pixels distribution and the charachteristic size of the galaxy

        float xmin=1e10;
        float xmax=-1e10;
        float ymin=1e10;
        float ymax=-1e10;
        float xm;
        float ym;
        float ref_size;
        float * max_mask = new float[nx*ny];
        int maxcount=0;
        for(long i; i<nx*ny; i++)max_mask[i]=0.0;
        for (int iy=1; iy<ny-1; iy++)
        for (int ix=1; ix<nx-1; ix++)
           {
             long index = ix+iy*nx;
             float IIImax=-0.5;
             for(int iiy=0; iiy<3; iiy++)
             for(int iix=0; iix<3; iix++)
               {
                   long index_aux = (ix-iix+1)+(iy-iiy+1)*nx;
                   IIImax = max(IIImax,III[index_aux]);
               }
             if(IIImax == III[index] && IIImax > 0.75)
               {
                   max_mask[index] = 1.0; 
                   cout << ix<< " " << iy << " " << max_mask[index] << endl;
                   maxcount++;
               }

             if(III[index] > 0.1)
               {
                   xmax = max(float(ix),xmax);
                   ymax = max(float(iy),xmax);
                   xmin = min(float(ix),xmin);
                   ymin = min(float(iy),xmin);
               }
               xm=xmax-xmin;
               ym=ymax-ymin; 
               ref_size = min(fabs(xm),fabs(ym));    
           }

// clean the maxima distribution
       
        int stencil = 3;
        for (int iy=stencil; iy<ny-stencil; iy++)
        for (int ix=stencil; ix<nx-stencil; ix++)
           {

             long index = ix+iy*nx;
             if (max_mask[index] == 1.0)
               {
                 for(int iiy=0; iiy<2*stencil+1; iiy++)
                 for(int iix=0; iix<2*stencil+1; iix++)
                   {
                      long index_aux = (ix-iix+stencil)+(iy-iiy+stencil)*nx;
                      if (max_mask[index_aux] == 1.0 && index != index_aux)
                        {
                           if(III[index] >= III[index_aux])max_mask[index_aux]=0;
                           if(III[index] < III[index_aux])max_mask[index]=0;
                        }
                   }
               }
           }
// create random (x,y,z) point distributions around each point:
// create clusters around III local maxima
// create diffuse component around all bright pixels

        //long nfinal = npergroup*maxcount;
        long nfinal = npergroup*nx*ny;
        xx = new float [nfinal]; 
        yy = new float [nfinal]; 
        zz = new float [nfinal]; 
        float r_dist;
        long * p_map = new long [nfinal];
        for(long iii=0; iii<nfinal; iii++)p_map[iii]=-1;
        long pcounter = 0;
        long p_i=0;
        cout << "SIZE " << ref_size << endl;
        for (int iy=0; iy<ny; iy++)
        for (int ix=0; ix<nx; ix++)
           {
             long index = ix+iy*nx;
             float xcenter = float(ix-nx/2);
             float ycenter = float(iy-ny/2);
             float radius_min = sqrt(xcenter*xcenter+ycenter*ycenter);
        
             //float beta=1.0;
             //float sigma_L = -(1.0/(beta*beta*ref_size))*radius_min*radius_min+ref_size;
             //if(sigma_L < 0.0)sigma_L=0.0;
             //////float radius = box_uniform(0,4*sigma_L);
             //////float radius = sigma_L;

             long IIIth=0;
             float r_0;
             float x_ref;
             if(III[index] > 0.0 && max_mask[index] != 1.0) 
               {
                 //IIIth = long(50*III[index]+3); 
                 IIIth = 20;
                 r_0 = ref_size*0.4;
                 x_ref = 0.5;
                 r_dist=1.0;
               } else if (max_mask[index] == 1.0) {
                 IIIth = long(npergroup*III[index]+10); 
                 r_0 = ref_size*0.8;
                 x_ref = 0.0;
                 r_dist=sigma_fixed;
               }

             float x_x = radius_min / r_0;
             if(x_x < x_ref) x_x = x_ref; 
             float sigma_L = x_x * exp(-x_x) * r_0;
             float radius = box_muller(0,sigma_L);
             float zcenter = radius;
             float zavg=0.0;
             int navg=0;

             if(IIIth > 0)
             {
               p_map[p_i] = pcounter;
               p_i++;
               long nrandom=IIIth;   
               for (int ii=0;ii<nrandom;ii++)
               {
                 xx[pcounter]=box_muller(xcenter,r_dist);
                 yy[pcounter]=box_muller(ycenter,r_dist);
                 zz[pcounter]=box_muller(zcenter,r_dist);
                 pcounter++;
               }
             }
           }
           p_map[p_i]=pcounter;

           for (long ii=0; ii<pcounter; ii++) 
           {  
                coordx[ii] = xx[ii]/norm;
                coordy[ii] = yy[ii]/norm;
                coordz[ii] = zz[ii]/norm;
           }


// write data
         pFile = fopen("test.dat", "w");
	 for (long ii=0; ii<pcounter; ii++)
             fprintf(pFile,"%f %f %f\n",xx[ii],
                                        yy[ii],
                                        zz[ii]);
         fclose(pFile);

         delete [] xx;
         delete [] yy;
         delete [] zz;
         return pcounter;
}



