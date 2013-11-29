# include "Galaxy.h"
#define MAX(a,b) ((a < b) ?  (b) : (a))
#define MIN(a,b) ((a > b) ?  (b) : (a))

float box_muller(float m, float s);

void CalculateColours (paramfile &params, string ComponentsName,long npart, 
                       float * cred, float * cgreen, float * cblue, float * ciii,
                       float * Red, float * Green, float * Blue, float * III, float * xcoord, 
                       float * ycoord, long nx, long ny)
{

	float xaux, yaux, xcol;
	long ii, jj, iaux;
        float x_rand_max = (float) RAND_MAX;
	float xcolaux;


	float brightness;
	float brightness_fact;
	brightness_fact = params.find<float>("Brightness"+ComponentsName,1.0);
	brightness = 1.0/(float)npart;
	brightness = brightness_fact*pow(brightness, 0.3333333f);
	///brightness = brightness_fact*sqrt(brightness);
	cout << "--> Coloring " << ComponentsName.c_str() << " with brightness " << brightness << endl;

        float pixeltotirific = params.find<float>("PixelToArcsec"+ComponentsName,-1.0);
        float nxarcsec = float(nx)*pixeltotirific;
        float nyarcsec = float(ny)*pixeltotirific;

        float iiimax=-1.0;
        float cextrema[6];
        cextrema[0]=1e20;
        cextrema[2]=1e20;
        cextrema[4]=1e20;
        cextrema[1]=-1e20;
        cextrema[3]=-1e20;
        cextrema[5]=-1e20;

	for (long particlei=0; particlei<npart; particlei++)
	{

           if(pixeltotirific == -1.0)
           {
	     xaux = (0.5*(xcoord[particlei]+1.0)); 
	     yaux = (0.5*(ycoord[particlei]+1.0)); 
	     ii = (int) (xaux*nx);
	     jj = (int) (yaux*ny);
           } else {
             xaux = xcoord[particlei]+nxarcsec*0.5;
             yaux = ycoord[particlei]+nyarcsec*0.5;
             ii = (int) ((xaux/nxarcsec)*nx);
             jj = (int) ((yaux/nxarcsec)*ny);
           }

	   if(ii >= nx || ii < 0 || jj >=ny || jj < 0)
	   {
//              xcol = ((float)rand())/x_rand_max;
	      xcolaux = box_muller(0, 0.25);
	      xcol = brightness*fabs(xcolaux);
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
	      ciii[particlei]   = brightness*III[iaux];
	   }
	
           iiimax = MAX(iiimax, ciii[particlei]);
           cextrema[0]=MIN(cextrema[0],cred[particlei]);
           cextrema[2]=MIN(cextrema[2],cgreen[particlei]);
           cextrema[4]=MIN(cextrema[4],cblue[particlei]);
           cextrema[1]=MAX(cextrema[1],cred[particlei]);
           cextrema[3]=MAX(cextrema[3],cgreen[particlei]);
           cextrema[5]=MAX(cextrema[5],cblue[particlei]);

	}
        cout << "Maximum Brightness = " << iiimax << endl;
        cout << "Red color range = " << cextrema[0] << " , " << cextrema[1] << endl;
        cout << "Green color range = " << cextrema[2] << " , " << cextrema[3] << endl;
        cout << "Blue color range = " << cextrema[4] << " , " << cextrema[5] << endl;
	//for (long particlei=0; particlei<npart; particlei++)ciii[particlei] /= iiimax;
}
