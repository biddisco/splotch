/*
 * Copyright (c) 2004-2006  
 *              Martin Reinecke (1), Klaus Dolag (1)
 *               (1) Max-Plank-Institute for Astrophysics
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>

#include "arr.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"

using namespace std;
using namespace RAYPP;

#define NO_SORTBYHSML
#define NO_SORTBYVAL
#define NO_REVERSESORT
#define NO_STOP_AFTER_N 5
#define NO_FULLABSROB
#define NO_MINHSMLPIXEL
#define NO_TIPSY
#define NO_PROJECTION_OFF
#define NO_PATCH_MACH
#define GEOMETRY_FILE
#define NO_STARS_ABSORB_RED
#define NO_PATCH_STARS
#define COLOR_VECTOR

const float8 Pi=3.14159265358979323846;

struct PARTICLE
  {
  float4 x,y,z,r,ro,I,T;
#ifdef COLOR_VECTOR
    float4 T2,T3;
#endif
  int type;
  };

struct zcmp
  {
  int operator()(const PARTICLE &p1, const PARTICLE &p2)
    {
#ifdef SORTBYVAL
#ifdef REVERSESORT
    return p1.T<p2.T;
#else
    return p1.T>p2.T;
#endif
#else
#ifdef SORTBYHSML
    return p1.r>p2.r;
#else
    return p1.z>p2.z;
#endif
#endif
    }
  };

class exptable
  {
  private:
    int nexp;
    float8 maxexptab;
    float8 expfac;
    float8 *exptab;

  public:
    exptable (int nexp_, float8 maxexp_)
      : nexp(nexp_), maxexptab(maxexp_), expfac(nexp/maxexptab),
        exptab(new float8[nexp])
      {
      for (int m=0; m<nexp; ++m)
        exptab[m]=exp(m/expfac);
      }

    float8 operator() (float8 arg) const
      {
      arg *= expfac;
      int iarg=int(arg);
      if (iarg<0) return 1.;
      if (iarg>=nexp-1) return 0.;
      float8 frac=arg-iarg;
      return (1.-frac)*exptab[iarg]+frac*exptab[iarg+1];
      }
  };

int gadget_find_block(bifstream &file,const string &label)
  {
  int i;
  int4 blocksize=0, blksize;
  char blocklabel[5]={"    "};

  file.clear();
  file.rewind();

  while(!file.eof() && blocksize == 0)
    {
    file >> blksize;
    if (blksize != 8)
      {
	cout << "wrong structure in GADGET file: " << blksize << endl;
      exit(1);
      }
    else
      {
      file.get(blocklabel,4);
      blocklabel[4] = 0;
      i=3;
      while (blocklabel[i]==32)
	{
	  blocklabel[i]=0;
	  i--;
	}
      cout << "blocklabel: <" << blocklabel << "> --> <" << label << ">" << endl;
      file >> blocksize;
      cout << "blocksize: " << blocksize << endl;
      file.skip(4);
      if (label!=blocklabel)
	{
        file.skip(blocksize);
        blocksize=0;
        }
      }
    }
  return(blocksize-8);
  }

int gadget_read_header(bifstream &file, int *npart,double *massarr,
  double &time,double &redshift, int *npartall)
  {
  int blocksize = gadget_find_block (file,"HEAD");
  planck_assert (blocksize>0, "Header block not found");
  file.skip(4);
  file.get(npart,6);
  file.get(massarr,6);
  file >> time >> redshift;
  file.skip(8);
  file.get(npartall,6);

  //  cout << " HEADER " << npart[0]<< " , " << npartall[0] << " Patching !!" << endl;  
  //  npartall[0]=34012224;
  //  npartall[1]=34012224;
  return blocksize;
  }

int main (int argc, const char ** argv)
  {
const float8 Pi=3.14159265358979323846;
  announce ("splotch");
  planck_assert (argc==2,"usage: splotch <parameter file>");
  paramfile params (argv[1]);

  cout << "initializing" << endl;
  int res = params.find<int>("resolution",200);
  int ycut0 = params.find<int>("ycut0",0);
  int ycut1 = params.find<int>("ycut1",res);
  string infilename = params.find<string>("infile");
  string outfile = params.find<string>("outfile");
  string filename;
  int numfiles = params.find<int>("numfiles",1);
  bool doswap = params.find<bool>("swap_endian",true);

  string lable_int = params.find<string>("lable_intensity","XXXX");
  string lable_col = params.find<string>("lable_color","XXXX");

  bool boostcolors = params.find<bool>("boostcolors",false);

  bool log_int = params.find<bool>("log_intensity",true);
  bool log_col = params.find<bool>("log_color",true);
  float minval_int = params.find<float>("min_int",1024);
  float maxval_int = params.find<float>("max_int",1024);
  float minval_col = params.find<float>("min_col",1024);
  float maxval_col = params.find<float>("max_col",1024);
  float hsml_fac = params.find<float>("hsml_fac",1.0);

  bool stars = params.find<bool>("stars",false);
  float star_hsml = params.find<float>("star_hsml",5.0);
  string stars_hsml_lable = params.find<string>("stars_hsml","XXXX");

  string stars_int = params.find<string>("stars_intensity","XXXX");
  string stars_col = params.find<string>("stars_color","XXXX");

  bool log_stars_int = params.find<bool>("log_stars_intensity",true);
  bool log_stars_col = params.find<bool>("log_stars_color",true);
  float minval_stars_int = params.find<float>("min_stars_int",1024);
  float maxval_stars_int = params.find<float>("max_stars_int",1024);
  float minval_stars_col = params.find<float>("min_stars_col",1024);
  float maxval_stars_col = params.find<float>("max_stars_col",1024);
  
  float8 brightness_stars = params.find<double>("brightness_stars",1.);
  float8 brightness_gas = params.find<double>("brightness_gas",1.);

  exptable xexp(10000,-20.);

  int4 npart[6],npartall[6];
  float8 massarr[6],time,redshift;

  if(numfiles>1) filename=infilename+"."+dataToString(0);
  else           filename=infilename;
  bifstream infile(filename.c_str(),doswap);
  if (!infile) {cout << "could not open input file! <" << filename << ">" << endl; exit(1);}
  gadget_read_header(infile,npart,massarr,time,redshift,npartall);
  if(!stars) npartall[4]=0;
#ifdef PATCH_STARS
  if(stars) npartall[4]=npartall[2]+npartall[3];
#endif
  cout <<  time << " " << redshift << endl;
  infile.close();

  int np=npartall[0];
  int nstar=npartall[4];
  cout << "Allocating space for " << np << " gas particles and " << nstar 
       << " star particles ..." << endl;  
  PARTICLE *p=new PARTICLE[np+nstar];
#ifdef GEOMETRY_FILE
  PARTICLE *p_orig=new PARTICLE[np+nstar];
#endif
  float4 mintmp=1e30, maxtmp=-1e30,minint=1e30, maxint=-1e30;
#ifdef COLOR_VECTOR
  float4 mintmp2=1e30, maxtmp2=-1e30,mintmp3=1e30, maxtmp3=-1e30;
#endif
  float4 mintmp_stars=1e30, maxtmp_stars=-1e30,minint_stars=1e30, maxint_stars=-1e30;
  int nskip,npstart=0,nstarstart=np;

  for(int f=0;f<numfiles;f++)
    {
     if(numfiles>1) filename=infilename+"."+dataToString(f);
     else           filename=infilename;
     bifstream infile(filename.c_str(),doswap);
     if (!infile) {cout << "could not open input file! <" << filename << ">" << endl; exit(1);}
     gadget_read_header(infile,npart,massarr,time,redshift,npartall);
     if(!stars) npart[4]=0;
     nskip=npart[1]+npart[2]+npart[3];
#ifdef PATCH_STARS
     if(stars) npart[4]=npart[2]+npart[3];
     nskip=npart[1];
#endif
     cout << "Reading " << npart[0] << " particles + " << npart[4] << " stars " << endl;
     cout << "Distributing at " << npstart << " and " << nstarstart << endl;
     gadget_find_block(infile,"POS");
     infile.skip(4);
     for (int m=0; m<npart[0]; ++m)
       infile >> p[npstart+m].x >> p[npstart+m].y >> p[npstart+m].z;
     if(stars)
       {
       cout << "Skipping " << nskip << " Positions ..." << endl;
       if(nskip>0) infile.skip(nskip*4*3);
       if(npart[4]>0)
         {
         for (int m=0; m<npart[4]; ++m)
           infile >> p[m+nstarstart].x >> p[m+nstarstart].y >> p[m+nstarstart].z;
         }
       }

     if(npart[0] > 0)
       {
	 for (int m=0; m<npart[0]; ++m)
	   p[m+npstart].type=0;

         cout << "Reading HSML (gas) ..." << endl;
	 gadget_find_block(infile,"HSML");
	 infile.skip(4);
	 for (int m=0; m<npart[0]; ++m)
	   {
	     infile >> p[m+npstart].r;
	     p[m+npstart].r *= hsml_fac;
	   }

	 cout << "Reading Colors (gas) ..." << endl;
	 if(gadget_find_block(infile,lable_col) > 0)
	   {
	     infile.skip(4);
	     for (int m=0; m<npart[0]; ++m)
	       {
		 infile >> p[m+npstart].T;
#ifdef COLOR_VECTOR
		 infile >> p[m+npstart].T2 >> p[m+npstart].T3;
#endif
		 if(log_col) p[m+npstart].T = log(p[m+npstart].T);
		 mintmp=min(mintmp,p[m+npstart].T);
		 maxtmp=max(maxtmp,p[m+npstart].T);
#ifdef COLOR_VECTOR
		 mintmp2=min(mintmp2,p[m+npstart].T2);
		 maxtmp2=max(maxtmp2,p[m+npstart].T2);
		 mintmp3=min(mintmp3,p[m+npstart].T3);
		 maxtmp3=max(maxtmp3,p[m+npstart].T3);
#endif
	       }
	   }
	 else
	   {
	     for (int m=0; m<npart[0]; ++m)
	       p[m+npstart].T=1.0;
	     mintmp=1;
	     maxtmp=1;
	   }

	 cout << "Reading Intensity (gas) ..." << endl;
	 if(gadget_find_block(infile,lable_int) > 0)
	   {
	     infile.skip(4);
	     for (int m=0; m<npart[0]; ++m)
	       {
		 infile >> p[m+npstart].I;
		 if(log_int) p[m+npstart].I = log(p[m+npstart].I);
		 minint=min(minint,p[m+npstart].I);
		 maxint=max(maxint,p[m+npstart].I);
	       }
	   }
	 else
	   {
	     for (int m=0; m<npart[0]; ++m)
	       p[m+npstart].I=1.0;
	     minint=1;
	     maxint=1;
	   }

       }
     if(stars)
       {
       if(npart[4]>0)
         {
	   for (int m=0; m<npart[4]; ++m)
             p[m+nstarstart].type=1;

	   cout << "Reading HSML (stars) ..." << endl;
	   if(gadget_find_block(infile,stars_hsml_lable) > 0)
	     {
	       infile.skip(4);
	       for (int m=0; m<npart[4]; ++m)
		 infile >> p[m+nstarstart].r;
	     }
	   else
	     for (int m=0; m<npart[4]; ++m)
	       p[m+nstarstart].r=star_hsml;

	   cout << "Reading Colors (stars) ..." << endl;
	   if(gadget_find_block(infile,stars_int) > 0)
	     {
	       infile.skip(4);
	       for (int m=0; m<npart[4]; ++m)
		 {
		   infile >> p[m+nstarstart].I;
		   if(log_int) p[m+nstarstart].I = log(p[m+nstarstart].I);
		   minint_stars=min(minint_stars,p[m+nstarstart].I);
		   maxint_stars=max(maxint_stars,p[m+nstarstart].I);
		 }
	     }
	   else
	     {
	       for (int m=0; m<npart[4]; ++m)
		 p[m+nstarstart].I=1.0;
	       minint_stars=1;
	       maxint_stars=1;
	     }

	   cout << "Reading Intensities (stars) ..." << endl;
	   if(gadget_find_block(infile,stars_col) > 0)
	     {
	       infile.skip(4);
	       for (int m=0; m<npart[4]; ++m)
		 {
		   infile >> p[m+nstarstart].T;
		   if(log_stars_col) p[m+nstarstart].T = log(p[m+nstarstart].T);
		   mintmp_stars=min(mintmp_stars,p[m+nstarstart].T);
		   maxtmp_stars=max(maxtmp_stars,p[m+nstarstart].T);
		 }
	     }
	   else
	     {
	       for (int m=0; m<npart[0]; ++m)
		 p[m+npstart].I=1.0;
	       mintmp_stars=1;
	       maxtmp_stars=1;
	     }

	 }
       }

     npstart+=npart[0];
     nstarstart+=npart[4];
    }

  if(minval_col!=maxval_col) 
    {
    mintmp=minval_col;
    maxtmp=maxval_col;
#ifdef COLOR_VECTOR
    mintmp2=minval_col;
    maxtmp2=maxval_col;
    mintmp3=minval_col;
    maxtmp3=maxval_col;
#endif
    }

  if(minval_int!=maxval_int) 
    {
    minint=minval_int;
    maxint=maxval_int;
    }

  if(minval_stars_col!=maxval_stars_col) 
    {
    mintmp_stars=minval_stars_col;
    maxtmp_stars=maxval_stars_col;
    }

  if(minval_stars_int!=maxval_stars_int) 
    {
    minint_stars=minval_stars_int;
    maxint_stars=maxval_stars_int;
    }

  cout << "Gas particles: " << endl;
  cout << "Color Range:     " << mintmp << " (min) , " << 
                                 maxtmp << " (max) " << endl;
  cout << "Intensity Range: " << minint << " (min) , " << 
                                 maxint << " (max) " << endl;
  cout << "Star particles: " << endl;
  cout << "Color Range:     " << mintmp_stars << " (min) , " << 
                                 maxtmp_stars << " (max) " << endl;
  cout << "Intensity Range: " << minint_stars << " (min) , " << 
                                 maxint_stars << " (max) " << endl;

  for (int m=0; m<np+nstar; ++m)
    {
    if(p[m].type==0)
      {
        if(mintmp != maxtmp) p[m].T = (p[m].T-mintmp)/(maxtmp-mintmp);
        if(minint != maxint) p[m].I = (p[m].I-minint)/(maxint-minint);
#ifdef COLOR_VECTOR
        if(mintmp2 != maxtmp2) p[m].T2 = (p[m].T2-mintmp2)/(maxtmp2-mintmp2);
        if(mintmp3 != maxtmp3) p[m].T3 = (p[m].T3-mintmp3)/(maxtmp3-mintmp3);
#endif
      }
    else
      {
        if(mintmp_stars != maxtmp_stars) p[m].T = (p[m].T-mintmp_stars)/(maxtmp_stars-mintmp_stars);
        if(minint_stars != maxint_stars) p[m].I = (p[m].I-minint_stars)/(maxint_stars-minint_stars);
      }
    }

  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*Pi/180.);

  COLOURMAP amap_gas,emap_gas,amap_stars,emap_stars;
  COLOURMAP *amap,*emap;

  FILE * pFile;
  string palname;
  int nColours;
  char foo[100];
  float dlut;
  float startlut;
  float rrr,ggg,bbb,rrr_old,ggg_old,bbb_old;

  for(int i=0;i<2;i++)
   {
     if (i==0) 
       {
         amap=&amap_gas;
         emap=&emap_gas;  
       }
     else
       {
         amap=&amap_stars;
         emap=&emap_stars;  
       }
     if(i==0)
       palname = params.find<string>("filegas_palette");
     else
       palname = params.find<string>("filestars_palette");

     pFile=fopen(palname.c_str(), "r");
     fscanf(pFile,"%s\n", foo ); // skip lut ide     == JASC-PAL
     fscanf(pFile,"%s\n", foo ); // skip lut version == 0100
     fscanf(pFile,"%d\n", &nColours ); 
     
     startlut = 0.0;
     dlut = 1.0/(float)(nColours-1);

     cout << "Loading " << nColours << " entries of color table"  << endl;             
     fscanf(pFile,"%f %f %f\n", &rrr_old, &ggg_old, &bbb_old );
         
     for (int i=1; i<nColours; i++)
       {
	 fscanf(pFile,"%f %f %f\n", &rrr, &ggg, &bbb );
	 amap->Add_Entry(new LINEAR_CMAP_ENTRY(startlut,startlut+dlut,
					       COLOUR(rrr_old/255,ggg_old/255,bbb_old/255),
					       COLOUR(rrr/255,ggg/255,bbb/255)));
	 rrr_old=rrr;
	 ggg_old=ggg;
	 bbb_old=bbb;  
	 startlut += dlut;
       }
     fclose(pFile);
     if(startlut < 1)
	 amap->Add_Entry(new LINEAR_CMAP_ENTRY(startlut,1,
					       COLOUR(rrr_old/255,ggg_old/255,bbb_old/255),
					       COLOUR(rrr/255,ggg/255,bbb/255)));
    }

  emap_gas=amap_gas;
  emap_stars=amap_stars;

  bool gas = params.find<bool>("gas",true);
  float zmaxval = params.find<float>("zmax",1.e23);
  float zminval = params.find<float>("zmin",0.0);
  float grayabsorb_gas = params.find<float>("gray_absorbtion_gas",0.2);
  float grayabsorb_stars = params.find<float>("gray_absorbtion_stars",0.2);

#ifdef GEOMETRY_FILE
  for (int m=0; m<np+nstar; ++m)
    {
       p_orig[m].x=p[m].x;
       p_orig[m].y=p[m].y;
       p_orig[m].z=p[m].z;
       p_orig[m].r=p[m].r;
       p_orig[m].ro=p[m].ro;
       p_orig[m].I=p[m].I;
       p_orig[m].T=p[m].T;
#ifdef COLOR_VECTOR
       p_orig[m].T2=p[m].T2;
       p_orig[m].T3=p[m].T3;
#endif
       p_orig[m].type=p[m].type;
    }
  float cam_x,cam_y,cam_z,lat_x,lat_y,lat_z,sky_x,sky_y,sky_z;
  string geometry_file = params.find<string>("geometry_file");
  string line;
  ifstream inp(geometry_file.c_str());
  int linecount=10000;
  int geometry_skip = params.find<int>("geometry_start",0);

  for(int i=0; i<geometry_skip; i++, linecount++)
    getline(inp, line);

  while (getline(inp, line))
    {
      for (int m=0; m<np+nstar; ++m)
        {
           p[m].x=p_orig[m].x;
           p[m].y=p_orig[m].y;
           p[m].z=p_orig[m].z;
           p[m].r=p_orig[m].r;
           p[m].ro=p_orig[m].ro;
           p[m].I=p_orig[m].I;
           p[m].T=p_orig[m].T;
#ifdef COLOR_VECTOR
           p[m].T2=p_orig[m].T2;
           p[m].T3=p_orig[m].T3;
#endif
           p[m].type=p_orig[m].type;
        }
    sscanf(line.c_str(),"%f %f %f %f %f %f %f %f %f",
                      &cam_x,&cam_y,&cam_z,
	              &lat_x,&lat_y,&lat_z,
                      &sky_x,&sky_y,&sky_z);
    VECTOR campos(cam_x,cam_y,cam_z);
    VECTOR lookat(lat_x,lat_y,lat_z);
    VECTOR sky(sky_x,sky_y,sky_z);
    cout << "Camera: (" << cam_x << "," << cam_y << "," << cam_z << ")" << endl;
    cout << "Lookat: (" << lat_x << "," << lat_y << "," << lat_z << ")" << endl;
    cout << "Sky:    (" << sky_x << "," << sky_y << "," << sky_z << ")" << endl;
#else
  VECTOR campos(params.find<double>("camera_x"),
                params.find<double>("camera_y"),
                params.find<double>("camera_z"));
  VECTOR lookat(params.find<double>("lookat_x"),
                params.find<double>("lookat_y"),
                params.find<double>("lookat_z"));
  VECTOR sky(params.find<double>("sky_x",0),
                params.find<double>("sky_y",1),
                params.find<double>("sky_z",0));
#endif
  cout << "transforming" << endl;
  sky.Normalize();
  VECTOR zaxis = (lookat-campos).Norm();
  VECTOR xaxis = Cross (sky,zaxis).Norm();
  VECTOR yaxis = Cross (zaxis,xaxis);
  TRANSFORM trans;
  trans.Make_General_Transform
    (TRANSMAT(xaxis.x,xaxis.y,xaxis.z,
              yaxis.x,yaxis.y,yaxis.z,
              zaxis.x,zaxis.y,zaxis.z,
              0,0,0));
  trans.Invert();
  TRANSFORM trans2;
  trans2.Make_Translation_Transform(-campos);
  trans2.Add_Transform(trans);
  trans=trans2;
  //  trans.Add_Transform(trans2);

  for (int m=0; m<np+nstar; ++m)
    {
    VECTOR v(p[m].x,p[m].y,p[m].z);
    v=trans.TransPoint(v);
    p[m].x=v.x; p[m].y=v.y; p[m].z=v.z;
    }
#ifdef PROJECTION_OFF
    float8 dist=sqrt((campos.x-lookat.x) * (campos.x-lookat.x)
                    +(campos.y-lookat.y) * (campos.y-lookat.y)
		    +(campos.z-lookat.z) * (campos.z-lookat.z));
    float8 xfac=1./(fovfct*dist);
    cout << "Field of fiew: " << 1./xfac*2. << endl;
#else
    float8 xfac;
#endif


  for (int m=0; m<np+nstar; ++m)
    {
#ifdef PROJECTION_OFF
    p[m].x = res*.5 * (p[m].x+fovfct*dist)*xfac;
    p[m].y = res*.5 * (p[m].y+fovfct*dist)*xfac;
#else
    xfac=1./(fovfct*p[m].z);
    p[m].x = res*.5 * (p[m].x+fovfct*p[m].z)*xfac;
    p[m].y = res*.5 * (p[m].y+fovfct*p[m].z)*xfac;
#endif
    p[m].ro = p[m].r;
    p[m].r = p[m].r *res*.5*xfac;
#ifdef MINHSMLPIXEL
    if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
      {
         p[m].r = 0.5;
         p[m].ro = p[m].r/res/0.5/xfac;
      }
#endif
    }

#ifdef SORTBYVAL
#ifdef REVERSESORT
  cout << "reverse sorting by value" << endl;
#else
  cout << "sorting by value" << endl;
#endif
#else
#ifdef SORTBYHSML
  cout << "sorting by hsml" << endl;
#else
  cout << "sorting by z" << endl;
#endif
#endif
  sort(&p[0], &p[np+nstar], zcmp());

#ifdef TIPSY
  cout << "plotting (adaptive smoothed tipsy like)" << endl;
#else
  cout << "plotting (raytracing)" << endl;
#endif
  if(boostcolors)
    cout << "Boosting colors when adding ..." << endl;

  float brightness,grayabsorb;

  float8 rfac=1.5;
  float8 powtmp = pow(float8(Pi/2),float8(1./3.));
  float8 sigma0=powtmp/sqrt(2*Pi);
  float8 i00,bfak=1./(2*sqrt(Pi)*powtmp);

  arr2<COLOUR> pic(res,res);
  for(int ix=0;ix<res;ix++)
	for(int iy=0;iy<res;iy++)
	  {
	    pic[ix][iy].r=0;
	    pic[ix][iy].g=0;
	    pic[ix][iy].b=0;
	  }

#ifdef STOP_AFTER_N
  int z_count[res][res];
#endif

  for (int m=0; m<np+nstar; ++m)
    {
    if ((m%100000)==0)
      cout << m << " of " << np+nstar << " Particles ..." << endl;

    if (p[m].z<=0) continue;
    if (p[m].z<=zminval) continue;
    if (p[m].z>=zmaxval) continue;
    //        cout << zminval << " " << p[m].z << " " << zmaxval << endl;
    if (!gas and p[m].type==0) continue;
    if (!stars and p[m].type==1) continue;
#ifdef STOP_AFTER_N
    int pxt=int(p[m].x),pyt=int(p[m].y);
    if (pxt < 0 || pxt >= res) continue;  
    if (pyt < 0 || pyt >= res) continue;  
    z_count[pxt][pyt]++;
    if (z_count[pxt][pyt] > STOP_AFTER_N) continue;
#endif

    if (p[m].type==1)
      {
         amap=&amap_stars;
         emap=&emap_stars;
         brightness=brightness_stars;       
         grayabsorb=grayabsorb_stars;       
      }
    else
      {
         amap=&amap_gas;
         emap=&emap_gas;
         brightness=brightness_gas;       
         grayabsorb=grayabsorb_gas;       
      }
    i00=brightness*bfak;

    float8 temp=p[m].T;
    temp=max(float8(0.0000001),min(float8(0.9999999),temp));
#ifdef COLOR_VECTOR
    float8 temp2=p[m].T2,temp3=p[m].T3;
    temp2=max(float8(0.0000001),min(float8(0.9999999),temp2));
    temp3=max(float8(0.0000001),min(float8(0.9999999),temp3));
#endif
    float8 intensity=p[m].I;
    intensity=max(float8(0.0000001),min(float8(0.9999999),intensity));
    COLOUR e;
#ifndef COLOR_VECTOR
    e=amap->Get_Colour(temp)*intensity;
#else
    if (p[m].type==1)
      e=amap->Get_Colour(temp)*intensity;
    else
      {
	e.r=temp*intensity;
	e.g=temp2*intensity;
	e.b=temp3*intensity;
      }
#endif
    float8 r=p[m].r;
    float8 radsq = rfac*rfac*r*r;
#ifndef TIPSY
    COLOUR a;
#ifdef FULLABSROB
    a=COLOUR(1,1,1);
#else
    a=e;
#endif
#ifdef STARS_ABSORB_RED
    if (p[m].type==1)  a.r=intensity;
#endif

    COLOUR q (e.r/(a.r+grayabsorb),e.g/(a.g+grayabsorb),e.b/(a.b+grayabsorb));

    float8 ro=p[m].ro;
    float8 sigma=ro*sigma0;
    float8 fac0=-0.5/(sigma*sigma);
    float8 i0=i00/ro;
    float8 prefac1 = ro/r*ro/r*fac0;
    float8 prefac2 = -0.5*i0;
#endif

    float8 posx=p[m].x, posy=p[m].y;
    int minx = max(int(posx-rfac*r),0), maxx = min(int(posx+rfac*r+1),res),
        miny = max(int(posy-rfac*r),ycut0), maxy = min(int(posy+rfac*r+1),ycut1);
    for (int y=miny; y<maxy; ++y)
      {
      float8 ysq=(y-posy)*(y-posy);
      for (int x=minx; x<maxx; ++x)
        {
        float8 dsq = (x-posx)*(x-posx) + ysq;
        if (dsq<radsq)
          {
#ifdef TIPSY
	    float8 uu=max(float8(0.),min(float8(1.0),sqrt(dsq/radsq)));
            float8 lfac=0.0;
            if(uu < 0.5) lfac=1.-6.*uu*uu*(1.-uu); 
            else         lfac=2*(1-uu)*(1-uu)*(1-uu);
	    pic[x][y].r = min(e.r * lfac + (1 - lfac) * pic[x][y].r,float8(1.0));
	    pic[x][y].g = min(e.g * lfac + (1 - lfac) * pic[x][y].g,float8(1.0));
	    pic[x][y].b = min(e.b * lfac + (1 - lfac) * pic[x][y].b,float8(1.0));
	    //	    	    pic[x][y].r = e.r;
	    //	    	    pic[x][y].g = e.g;
	    //	    	    pic[x][y].b = e.b;
#else
          float8 fac = prefac2*xexp(prefac1*dsq);
	  if(boostcolors)
	    {
	      pic[x][y].r = e.r+(pic[x][y].r-e.r)*xexp(fac*q.r);
	      pic[x][y].g = e.g+(pic[x][y].g-e.g)*xexp(fac*q.g);
	      pic[x][y].b = e.b+(pic[x][y].b-e.b)*xexp(fac*q.b);
	    }
	  else
	    {
	      pic[x][y].r = q.r+(pic[x][y].r-q.r)*xexp(fac*a.r);
	      pic[x][y].g = q.g+(pic[x][y].g-q.g)*xexp(fac*a.g);
	      pic[x][y].b = q.b+(pic[x][y].b-q.b)*xexp(fac*a.b);
	    }
#endif
          }
        }
      }
    }

  bool colbar = params.find<bool>("colourbar",false);
  if (colbar)
    {
      cout << "adding colour bar ..." << endl;
      for (int x=0; x<res; x++)
        {
           float8 temp=x/float8(res);
           COLOUR e=amap_gas.Get_Colour(temp);
           for (int y=0; y<10; y++)
	     {
	        pic[x][res-1-y].r = e.r;
    	        pic[x][res-1-y].g = e.g;
	        pic[x][res-1-y].b = e.b;
	     }
        }
      for (int x=0; x<res; x++)
        {
           float8 temp=x/float8(res);
           COLOUR e=amap_stars.Get_Colour(temp);
           for (int y=10; y<20; y++)
	     {
	        pic[x][res-1-y].r = e.r;
    	        pic[x][res-1-y].g = e.g;
	        pic[x][res-1-y].b = e.b;
	     }
        }
    }

  cout << "writing" << endl;
  int yres=ycut1-ycut0;
  const byte header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    res%256, res/256, yres%256, yres/256, 24, 32 };

#ifdef GEOMETRY_FILE
  char frame_name[300];
  sprintf(frame_name,"%s_%5i.tga",outfile.c_str(),linecount);
  bofstream file(frame_name,file_is_natural);
  linecount++;
#else
  bofstream file(outfile.c_str(),file_is_natural);
#endif
  for (int m=0; m<18; ++m) file << header[m];
  for (int y=ycut0; y<ycut1; ++y)
    {
    for (int x=0; x<res; ++x)
      {
      pic[x][y].b=min(float8(1.), max(float8(0.), float8(pic[x][y].b)));
      byte pix = min(byte(255),byte(256*pic[x][y].b));
      file << pix;
      pic[x][y].g=min(float8(1.), max(float8(0.), float8(pic[x][y].g)));
      pix = min(byte(255),byte(256*pic[x][y].g));
      file << pix;
      pic[x][y].r=min(float8(1.), max(float8(0.), float8(pic[x][y].r)));
      pix = min(byte(255),byte(256*pic[x][y].r));
      file << pix;
      }
    }

#ifdef GEOMETRY_FILE
    }
#endif
  }
