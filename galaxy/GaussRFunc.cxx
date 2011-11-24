# include "Galaxy.h"

float box_muller(float m, float s);

long GaussRFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, 
		 float * coordx, float * coordy, float * coordz,
                 float xmax, float ymax, float zmax, float * zdeep, float * III, long nx, long ny) 
{
	FILE * pFile;
	float mean = 0.0;
        float sigma[3];
	long npergroup = 8;
        long freeparts = 0;
	long norig = 0;
	long ns=0;
	long rx,ry;
	float * resolution;
	long nnnn = number_of_points;
	float compression;
        float gsigmax;
        float gsigmay;
	float sigma_aux;

	
        srand(time(NULL));

        int grouping = params.find<int>("Grouping"+ComponentName,0);
	if(grouping == 0)
	{
        sigma[0] = params.find<float>("Sigmax"+ComponentName,0);
        sigma[1] = params.find<float>("Sigmay"+ComponentName,0);
        sigma[2] = params.find<float>("Sigmaz"+ComponentName,0);

        printf("========================\n");
        printf("MEAN  = %f\n", mean);
        printf("SIGMA X = %f\n", sigma[0]);
        printf("SIGMA Y = %f\n", sigma[1]);
        printf("SIGMA Z = %f\n", sigma[2]);
        printf("Number of Points = %d\n", number_of_points);
        printf("========================\n");
	}

	float radius = 0.0;

	if(grouping == 1)
	{
           sigma[2] = params.find<float>("Sigmaz"+ComponentName,0);
           npergroup = params.find<long>("NperGroup"+ComponentName,0);
           rx = params.find<long>("Scaledxres"+ComponentName,1);
           ry = params.find<long>("Scaledyres"+ComponentName,1);
           compression = params.find<float>("CompressionFactor"+ComponentName,1.0);

	   //rx++;
	   //ry++;

	   resolution = new float [rx*ry];
	   long countin = 0;
	   long countout = 0;

	   for(long i=0; i<rx*ry; i++)
	      resolution[i] = 0;

	   gsigmax = 3.0/(float)rx;
	   gsigmay = 3.0/(float)ry;
           for (long i=0; i<number_of_points; i++)
           {

	     //coordx[i] = box_muller(coordx[i], gsigmax);
	     //coordy[i] = box_muller(coordy[i], gsigmay);
             int ix = (int) round((0.5*(coordx[i]+1.0))*(rx));
             int iy = (int) round((0.5*(coordy[i]+1.0))*(ry));
             int irx = (int) round((0.5*(coordx[i]+1.0))*nx);
             int iry = (int) round((0.5*(coordy[i]+1.0))*ny);
/*
             int ix =  (int) round((coordx[i]+0.5)*rx);
             int iy =  (int) round((coordy[i]+0.5)*ry);
             int irx = (int) round((coordx[i]+0.5)*nx);
             int iry = (int) round((coordy[i]+0.5)*ny);
*/
	     

	     coordz[i]=100000.0;

	     long index = ix+iy*rx;
	     long rindex = irx+iry*nx;

	     if(iy <= ry || iy >= 0)
	     {

	       if(resolution[index] == 0.0)
	       {
		  float compp; 
		  if(III[rindex] >= 0.8)
		    {
			compp = 6.5*III[rindex];
		    }else{
			compp = 1.0;
		    }
		  sigma_aux = sigma[2]/compp;
/* KLAUS
                  sigma_aux=f(signal[rindex]); where f is an appropriate function
*/
		  coordz[i] = box_muller(mean, sigma_aux);
                  resolution[index] = 1.0;
		  countin++;
	       } else {
//	          coordz[i] = resolution[index];
		  coordx[i] = 0.0;
		  coordy[i] = 0.0;
		  coordz[i] = 0.0;
                  countout++;

	       }
	     }
           }

	   float * xaux = new float [countin];
	   float * yaux = new float [countin];
	   float * zaux = new float [countin];
	   long iaux = 0;

	   for (long i=0; i<number_of_points; i++)
              {
		if(coordx[i] == 0.0 && coordy[i] == 0.0 && coordz[i] == 0.0)continue;
		xaux[iaux] = coordx[i];
		yaux[iaux] = coordy[i];
		zaux[iaux] = coordz[i];
		iaux++;
		coordx[i] = 0.0;
		coordy[i] = 0.0;
		coordz[i] = 0.0;
	      }


	   for (long i=0; i<countin; i++)
              {
		coordx[i] = xaux[i];
		coordy[i] = yaux[i];
		coordz[i] = zaux[i];
	      }
	   delete [] xaux;
	   delete [] yaux;
	   delete [] zaux;

	   number_of_points = countin;



	} else {


// x coord

        for (long i=0; i<number_of_points; i++)
        {
             coordx[i] = box_muller(mean, sigma[0]);
        }


// y coord

        for (long i=0; i<number_of_points; i++)
        {
             coordy[i] = box_muller(mean, sigma[1]);
        }


// z coord


        for (long i=0; i<number_of_points; i++)
        {
	     
             coordz[i] = box_muller(mean, sigma[2]);
        }

	}


// grouping

	long abscounter=number_of_points;
	float gsigmaz = sigma[2]/compression;
//	gsigma = 0.1/(float)rx;
	long ii = 0;

	
	if(grouping == 1)
	{

	   while (abscounter < ntot-npergroup && ii < number_of_points)
	   {
	      //if(coordx[ii] != 0.0 && coordy[ii] != 0.0 && coordz[ii] != 0.0)
	      //{
	      for(long jj=0; jj<npergroup; jj++)
	      {
                 coordx[abscounter] = box_muller(coordx[ii], gsigmax);
                 coordy[abscounter] = box_muller(coordy[ii], gsigmay);
/* KLAUS
                 sigma_aux=f(signal[rindex]); where f is an appropriate function
                 gsigmaz = sigma_aux/compression; 
*/
                 coordz[abscounter] = box_muller(coordz[ii], gsigmaz);
		 abscounter++;
	      }
	      //}	      
              ii++;

	   }

	   nnnn = abscounter-1;
	}
		

        return nnnn;

}
