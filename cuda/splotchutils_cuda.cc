#include "cuda/splotchutils_cuda.h"

using namespace std;

/////////////////////////////CUDA CODE///////////////////////////////////
#ifdef CUDA
void render_as_thread1 (const vector<particle_sim> &p, arr2<COLOUR> &pic, 
      bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*pi);
      const float64 bfak=1./(2*sqrt(pi)*powtmp);
	  exptable xexp(-20.);

      int xres = pic.size1(), yres=pic.size2();
      pic.fill(COLOUR(0,0,0));

#ifdef VS
	  work_distributor wd (xres,yres,xres,yres);
#else
      work_distributor wd (xres,yres,200,200);
#endif //ifdef VS

    #pragma omp parallel
{
      int chunk;
    #pragma omp for schedule(dynamic,1)
      for (chunk=0; chunk<wd.nchunks(); ++chunk)
        {
        int x0, x1, y0, y1;
        wd.chunk_info(chunk,x0,x1,y0,y1);
        arr2<COLOUR> lpic(x1-x0,y1-y0);
        lpic.fill(COLOUR(0,0,0));
        int x0s=x0, y0s=y0;
        x1-=x0; x0=0; y1-=y0; y0=0;


        for (unsigned int m=0; m<p.size(); ++m)
	if(p[m].active==1)
          {
          float64 r=p[m].r;
          float64 posx=p[m].x, posy=p[m].y;
          posx-=x0s; posy-=y0s;
          float64 rfacr=rfac*r;

		  //in one chunk this culling is not necessary as it was done in coloring
          int minx=int(posx-rfacr+1);
          if (minx>=x1) continue;
          minx=max(minx,x0);
          int maxx=int(posx+rfacr+1);
          if (maxx<=x0) continue;
          maxx=min(maxx,x1);
          if (minx>=maxx) continue;
          int miny=int(posy-rfacr+1);
          if (miny>=y1) continue;
          miny=max(miny,y0);
          int maxy=int(posy+rfacr+1);
          if (maxy<=y0) continue;
          maxy=min(maxy,y1);
          if (miny>=maxy) continue;

	  COLOUR8 a=p[m].e, e, q;
          if (!a_eq_e)
            {
            e=p[m].e;
            q=COLOUR8(e.r/(a.r+grayabsorb),e.g/(a.g+grayabsorb),e.b/(a.b+grayabsorb));
            }

          float64 radsq = rfacr*rfacr;
          float64 prefac1 = -0.5/(r*r*sigma0*sigma0);
          float64 prefac2 = -0.5*bfak/p[m].ro;
          for (int x=minx; x<maxx; ++x)
            {
            float64 xsq=(x-posx)*(x-posx);
            for (int y=miny; y<maxy; ++y)
              {
              float64 dsq = (y-posy)*(y-posy) + xsq;
              if (dsq<radsq)
                {
                float64 fac = prefac2*xexp(prefac1*dsq);
                if (a_eq_e)
                  {
                  lpic[x][y].r += (fac*a.r);
                  lpic[x][y].g += (fac*a.g);
                  lpic[x][y].b += (fac*a.b);
                  }
                else
                  {
                  lpic[x][y].r = q.r+(lpic[x][y].r-q.r)*xexp(fac*a.r);
                  lpic[x][y].g = q.g+(lpic[x][y].g-q.g)*xexp(fac*a.g);
                  lpic[x][y].b = q.b+(lpic[x][y].b-q.b)*xexp(fac*a.b);
                  }//if a_eq_e
                }// if dsq<radsq
              }//y
            }//x
          }//for particle[m]
        for(int ix=0;ix<x1;ix++)
          for(int iy=0;iy<y1;iy++)
            pic[ix+x0s][iy+y0s]=lpic[ix][iy];
        }//for this chunk
}//#pragma omp parallel
}

#ifdef HOST_THREAD_RENDER
//DWORD WINAPI render_thread (const vector<particle_splotch> &p, int start, int end, arr2<COLOUR> &pic, 
//      bool a_eq_e,double grayabsorb)
DWORD WINAPI render_thread (param_render_thread *param)
      {
		cu_particle_splotch *p =param->p;
//		cu_color pic[][800] =param->pic;
	    bool a_eq_e =param->a_eq_e;
		double grayabsorb =param->grayabsorb;

      const float64 rfac=1.5;
      const float64 powtmp = pow(pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*pi);
      const float64 bfak=1./(2*sqrt(pi)*powtmp);
	  exptable xexp(-20.);

      int xres = 800, yres=800;
      memset(pic, 0, sizeof(cu_color) *800 *800);

	  work_distributor wd (xres,yres,xres,yres);

    #pragma omp parallel
{
      int chunk;
    #pragma omp for schedule(dynamic,1)
      for (chunk=0; chunk<wd.nchunks(); ++chunk)
        {
        int x0, x1, y0, y1;
        wd.chunk_info(chunk,x0,x1,y0,y1);
        arr2<COLOUR> lpic(x1-x0,y1-y0);
        lpic.fill(COLOUR(0,0,0));
        int x0s=x0, y0s=y0;
        x1-=x0; x0=0; y1-=y0; y0=0;


        for (unsigned int m=param->start; m<param->end; ++m)
          {
          float64 r=p[m].r;
          float64 posx=p[m].x, posy=p[m].y;
          posx-=x0s; posy-=y0s;
          float64 rfacr=rfac*r;

		  //in one chunk this culling is not necessary as it was done in coloring
          int minx=int(posx-rfacr+1);
          if (minx>=x1) continue;
          minx=max(minx,x0);
          int maxx=int(posx+rfacr+1);
          if (maxx<=x0) continue;
          maxx=min(maxx,x1);
          if (minx>=maxx) continue;
          int miny=int(posy-rfacr+1);
          if (miny>=y1) continue;
          miny=max(miny,y0);
          int maxy=int(posy+rfacr+1);
          if (maxy<=y0) continue;
          maxy=min(maxy,y1);
          if (miny>=maxy) continue;

	  cu_color a=p[m].e, e, q;
          if (!a_eq_e)
            {
            e=p[m].e;
			q.r=e.r/(a.r+grayabsorb);
			q.g=e.g/(a.g+grayabsorb);
			q.b=e.b/(a.b+grayabsorb);
            }

          float64 radsq = rfacr*rfacr;
          float64 prefac1 = -0.5/(r*r*sigma0*sigma0);
          float64 prefac2 = -0.5*bfak/p[m].ro;
          for (int x=minx; x<maxx; ++x)
            {
            float64 xsq=(x-posx)*(x-posx);
            for (int y=miny; y<maxy; ++y)
              {
              float64 dsq = (y-posy)*(y-posy) + xsq;
              if (dsq<radsq)
                {
                float64 fac = prefac2*xexp(prefac1*dsq);
                if (a_eq_e)
                  {
                  pic[x][y].r += (fac*a.r);
                  pic[x][y].g += (fac*a.g);
                  pic[x][y].b += (fac*a.b);
                  }
                else
                  {
                  pic[x][y].r = q.r+(lpic[x][y].r-q.r)*xexp(fac*a.r);
                  pic[x][y].g = q.g+(lpic[x][y].g-q.g)*xexp(fac*a.g);
                  pic[x][y].b = q.b+(lpic[x][y].b-q.b)*xexp(fac*a.b);
                  }//if a_eq_e
                }// if dsq<radsq
              }//y
            }//x
          }//for particle[m]
        }//for this chunk
}//#pragma omp parallel
		
		//post-process is deleted
		return (param->end - param->start);
}
#endif //ifdef HOST_THREAD_RENDER
#endif //ifdef CUDA
#ifdef CUDA
extern "C" 
void getCuTransformParams(cu_param_transform &para_trans,
paramfile &params, double c[3], double l[3], double s[3])
{
//	cout<<endl<<"\nRetrieve parameters for device transformation\n"<<endl;

	vec3	campos(c[0], c[1], c[2]),
			lookat(l[0], l[1], l[2]), 
			sky(s[0], s[1], s[2]);
	
	int res = params.find<int>("resolution",200);
	double fov = params.find<double>("fov",45); //in degrees
	double fovfct = tan(fov*0.5*degr2rad);
	float64 xfac=0.0, dist=0.0;

	sky.Normalize();
	vec3 zaxis = (lookat-campos).Norm();
	vec3 xaxis = crossprod (sky,zaxis).Norm();
	vec3 yaxis = crossprod (zaxis,xaxis);
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
	bool projection = params.find<bool>("projection",true);

	if(!projection)
	{
	  dist= (campos-lookat).Length();
	  xfac=1./(fovfct*dist);
	  cout << " Field of fiew: " << 1./xfac*2. << endl;
	}

	bool minhsmlpixel = params.find<bool>("minhsmlpixel",false);

	//retrieve the parameters for tansformation
	for (int i=0; i<12; i++)
		para_trans.p[i] =trans.Matrix().p[i];
	para_trans.projection	=projection;
	para_trans.res			=res;
	para_trans.fovfct		=fovfct;
	para_trans.dist			=dist;
	para_trans.xfac			=xfac;
	para_trans.minhsmlpixel =minhsmlpixel;
}
#endif	//CUDA
