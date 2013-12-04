# include "Galaxy.h"
#include"cxxsupport/arr.h"
#include"kernel/colourmap.h"
#include"cxxsupport/vec3.h"

/*

This function allows to warp an initial "regular" point distribution 
according to the Tirific model.
Conversion factor from original points units to arcsec (Tirific unit) is necessary
Final points coordinates are in arcsec

*/

// this is used to create spheroids of points: OPTION 1
long TirificWarp (paramfile &params, string ComponentName, long number_of_points, 
		 float * coordx, float * coordy, float * coordz) 
{

// in order to rotate, we need to go through the TIRIFIC stuff
// Now physical dimension comes into play: this is set by sigma
// if sigma is in arcsec then parsectotirific = 1.0
// if sigma is in kpc then parsectotirific must be properly set
//
// Read Tirific file (only radius Inclination and Position Angle actually needed)

  float pi=3.141592654;
  float parsectotirific = params.find<float>("InternToArcsec"+ComponentName,1.0);
  COLOURMAP model;

  ifstream infile (params.find<string>("TirificModel"+ComponentName).c_str());
  planck_assert (infile,"could not open TIRIFIC file  <" + params.find<string>("TirificModel"+ComponentName) + ">");
  string dummy;
  int nModel;
  infile >> nModel;
  infile >> dummy >> dummy >> dummy >> dummy >> dummy;
  cout << "      loading " << nModel << " entries of tirific model table " << endl;
  float rrr,vvv,zzz,iii,ppp;
  for (int i=0; i<nModel; i++)
    {
      infile >> rrr >> vvv >> zzz >> iii >> ppp;
      model.addVal(rrr,COLOUR(zzz,iii,ppp));
    }
   infile.close();

// Rotate particles

  long icount = 0;
  float r0 = 0.0;

  float xmax=-1e20;
  float xmin=1e20;
  float xcoord[3];
  long index;
  long ntrial = number_of_points;
  for (long i=0; i<ntrial; i++)    // loop over all possible particles
    {
      index = i;
      xcoord[0] = coordx[index];
      xcoord[1] = coordy[index];
      xcoord[2] = coordz[index];
      r0 = sqrt(xcoord[0]*xcoord[0] + xcoord[1]*xcoord[1]);

      double dphi = asin(double(xcoord[1]/r0));
      if(xcoord[1] >= 0.0 && xcoord[0] < 0.0)dphi = pi - dphi;
      if(xcoord[1] < 0.0 && xcoord[0] < 0.0)dphi = pi - dphi;
      if(xcoord[1] < 0.0 && xcoord[0] >= 0.0)dphi = 2.0*pi + dphi;
      float phi = float(dphi);

      float r1m = r0;     // radius associated to the particle
// Interpolate the tirific model to this radial distance
      COLOUR ring = model.getVal(r1m * parsectotirific);
      float thick = ring.r;
      float pa = (ring.b + 180) / 180 * M_PI;
      float inc = ring.g / 180 * M_PI;

      // Calculate the normalized normal vector of the ring from the tirific model
      vec3 nn(0,0,1);
      vec3 m1_1(1,0,0),m1_2(0,cos(inc),-sin(inc)),m1_3(0,sin(inc),cos(inc));
      vec3 m2_1(cos(pa),-sin(pa),0),m2_2(sin(pa),cos(pa),0),m2_3(0,0,1);
      vec3 n1(dotprod(m1_1,nn),dotprod(m1_2,nn),dotprod(m1_3,nn));
      vec3 n(dotprod(m2_1,n1),dotprod(m2_2,n1),dotprod(m2_3,n1));
      n.Normalize();

      // Find the radius vector within the x/y plane
      float x0;
      float y0;
      if (inc == 0.0)
      {
         x0=xcoord[0];
         y0=xcoord[1];
      } else {
         x0 = sqrt(r1m*r1m/(1+(n.x/n.y)*(n.x/n.y)));
         y0 = -(n.x/n.y)*x0;
      }
      vec3 rr(x0,y0,0);

      // Rotate the radius vector arround the normal vector by phi
      float sp = sin(phi);
      float cp = cos(phi);
      vec3 m3_1(n.x*n.x*(1-cp)+cp    ,n.x*n.y*(1-cp)-n.z*sp,n.x*n.z*(1-cp)+n.y*sp);
      vec3 m3_2(n.y*n.x*(1-cp)+n.z*sp,n.y*n.y*(1-cp)+cp    ,n.y*n.z*(1-cp)-n.x*sp);
      vec3 m3_3(n.z*n.x*(1-cp)-n.y*sp,n.z*n.y*(1-cp)+n.x*sp,n.z*n.z*(1-cp)+cp);
      vec3 xx(dotprod(m3_1,rr),dotprod(m3_2,rr),dotprod(m3_3,rr));

      float height = xcoord[2];
      vec3 xxx = xx + n * height;

      coordx[index] = xxx.x*parsectotirific;
      coordy[index] = xxx.y*parsectotirific;
      coordz[index] = xxx.z*parsectotirific;

    }
  return number_of_points;

}
