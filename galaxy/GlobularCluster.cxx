# include "Galaxy.h"

float box_muller(float m, float s);

long GlobularCluster (paramfile &params, string ComponentName, long number_of_points, long npergroup, 
		      float * coordx, float * coordy, float * coordz)
{
	float mean = 0.0;
        float sigma[3];
	long nnnn = number_of_points;
        float gsigmax;
	
        srand(time(NULL));

        gsigmax = params.find<float>("Gsigma",0);
        sigma[0] = params.find<float>("Sigmax"+ComponentName,0);
        sigma[1] = params.find<float>("Sigmay"+ComponentName,0);
        sigma[2] = params.find<float>("Sigmaz"+ComponentName,0);

        printf("========================\n");
        printf("MEAN  = %f\n", mean);
        printf("SIGMA X = %f\n", sigma[0]);
        printf("SIGMA Y = %f\n", sigma[1]);
        printf("SIGMA Z = %f\n", sigma[2]);
        printf("SIGMA G = %f\n", gsigmax);
        printf("Number of Points = %d\n", number_of_points);
        printf("========================\n");

	if(number_of_points > 0)
	{
// x coord

        for (long i=0; i<number_of_points; i++)
             coordx[i] = box_muller(mean, sigma[0]);

// y coord

        for (long i=0; i<number_of_points; i++)
             coordy[i] = box_muller(mean, sigma[1]);

// z coord

        for (long i=0; i<number_of_points; i++)
             coordz[i] = box_muller(mean, sigma[2]);

// grouping
/*
	long nincluster;
	long abscounter=number_of_points;
	for (long i=0; i<number_of_points; i++)
	    {

		nincluster = (long)((float)rand()/((float)RAND_MAX)*1000)+100;
		for (long j=0; j<nincluster; j++)
		    {

                 	coordx[abscounter] = box_muller(coordx[i], gsigmax);
                 	coordy[abscounter] = box_muller(coordy[i], gsigmax);
                 	coordz[abscounter] = box_muller(coordz[i], gsigmax);
		 	abscounter++;
	      	    }
	    }

	nnnn = abscounter-1;
*/
	}
        return nnnn;

}
