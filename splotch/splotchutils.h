#ifndef SPLOTCHUTILS_H
#define SPLOTCHUTILS_H

#ifdef VS
#include "cxxsupport/mpi_support.h"
#else
#include "mpi_support.h"
#endif
#include "kernel/transform.h"

#ifdef CUDA
	#include "cuda/splotch_cuda.h"

	//just for now that it's still needed for cuda codes 2 Dec 2009
	struct particle_splotch
	{
		float32 x,y,r,ro;
		COLOUR	a, e;

		particle_splotch (float32 x_, float32 y_, float32 r_, float32 ro_, const COLOUR &a_, const COLOUR &e_)
		: x(x_), y(y_), r(r_), ro(ro_), a(a_), e(e_) {}
		particle_splotch () {}
	};
	
	//#define CUDA_TEST_FRAGMENT
	
	#ifdef CUDA_TEST_FRAGMENT
		//the fragment buffer at host
		unsigned long	posFragBufH=0;
		cu_fragment_AeqE	fragBuf[30*100000*12];//for current test only! 10-29
		cu_fragment_AeqE	fragBufWrittenByHost[30*100000*12];//for current test only! 10-29
	#endif //ifdef CUDA_TEST_FRAGMENT

	//a array used for debug only
	#ifdef CUDA_TEST_COLORMAP
		int		size_inout_buffer=0;	//for debug only.
		float	inout_buffer[500000][5];//for debug only. ptype,input,r,g,b, 
	#endif //if CUDA_TEST_COLORMAP

	//#define CUDA_TEST_EXP
	#ifdef CUDA_TEST_EXP
		int		size_inout_buffer=0;	//for debug only.
		float	inout_buffer[100000][2];//for debug only. in, out
	#endif //if CUDA_TEST_EXP

	//#define HOST_THREAD_RENDER
	#ifdef HOST_THREAD_RENDER
	struct particle_splotch;
	struct	param_render_thread
	{
		cu_particle_splotch	*p;
		int start, end;
	    bool a_eq_e;
		double grayabsorb;
		cu_color pic[][800];
//		void	**pic;
	};
	#endif //ifdef HOST_THREAD_RENDER
#endif //if CUDA


using namespace std;
using namespace RAYPP;

class COLOUR8
  {
  public:
    float64 r, g, b;

    COLOUR8 () {}
    COLOUR8 (float64 rv, float64 gv, float64 bv)
      : r (rv), g (gv), b (bv) {}
    COLOUR8 (const COLOUR &orig)
      : r (orig.r), g (orig.g), b (orig.b) {}
    operator COLOUR()
      { return COLOUR (r,g,b); }

    COLOUR8 &operator=(const COLOUR &orig)
      {
      r=orig.r; g=orig.g; b=orig.b;
      return *this;
      }
  };

struct particle_sim
  {
  float32 x,y,z,r,ro,I,C1,C2,C3;
  int type,active;
  COLOUR e;

#ifdef INTERPOLATE
  unsigned int id;
#ifdef HIGH_ORDER_INTERPOLATION
  float32 vx,vy,vz;
#endif
#endif
  particle_sim (float32 x_, float32 y_, float32 z_, float32 r_, float32 ro_, 
                float32 I_, float32 C1_, float32 C2_, float32 C3_, int type_,
                int active_, const COLOUR &e_
#ifdef INTERPOLATE
                , unsigned int id_
#ifdef HIGH_ORDER_INTERPOLATION
                , float32 vx_, float32 vy_, float32 vz_
#endif
#endif
                ): x(x_), y(y_), z(z_), r(r_), ro(ro_), I(I_), C1(C1_), C2(C2_), C3(C3_), type(type_), active(active_), e(e_)
#ifdef INTERPOLATE
                , id(id_)
#ifdef HIGH_ORDER_INTERPOLATION
                , vx(vx_), vy(vy_), vz(vz_)
#endif
#endif
                 {}
  particle_sim () {}

  };


struct zcmp
  {
  int operator()(const particle_sim &p1, const particle_sim &p2)
    {
    return p1.z>p2.z;
    }
  };

struct vcmp1
  {
  int operator()(const particle_sim &p1, const particle_sim &p2)
    {
    return p1.C1>p2.C1;
    }
  };

struct vcmp2
  {
  int operator()(const particle_sim &p1, const particle_sim &p2)
    {
    return p1.C1<p2.C1;
    }
  };

struct hcmp
  {
  int operator()(const particle_sim &p1, const particle_sim &p2)
    {
    return p1.r>p2.r;
    }
  };

template<typename T> void get_minmax (T &minv, T &maxv, T val)
  {
  minv=min(minv,val);
  maxv=max(maxv,val);
  }

double my_asinh (double val)
  {
  return log(val+sqrt(1.+val*val));
  }

template<typename T> void my_normalize (T minv, T maxv, T &val)
  {
  if (minv!=maxv) val =  (val-minv)/(maxv-minv);
  }

template<typename T> void clamp (T minv, T maxv, T &val)
  {
  val = min(maxv, max(minv, val));
  }


class exptable
  {
  private:
    float64 expfac;
    arr<float64> tab1, tab2;
    enum {
      nbits=10,
      dim1=1<<nbits,
      mask1=dim1-1,
      dim2=1<<nbits<<nbits,
      mask2=dim2-1,
      mask3=~mask2
      };

  public:
    exptable (float64 maxexp)
      : expfac(dim2/maxexp), tab1(dim1), tab2(dim1)
      {
      for (int m=0; m<dim1; ++m)
        {
        tab1[m]=exp(m*dim1/expfac);
        tab2[m]=exp(m/expfac);
		}
      }

    float64 operator() (float64 arg) const
      {
      int iarg= (int)(arg*expfac);
      if (iarg&mask3)
        return (iarg<0) ? 1. : 0.;
      return tab1[iarg>>nbits]*tab2[iarg&mask1];
      }
  };

class work_distributor
  {
  private:
    int sx, sy, tx, ty;

  public:
    work_distributor (int sx_, int sy_, int tx_, int ty_)
      : sx(sx_), sy(sy_), tx(tx_), ty(ty_) {}

    int nchunks() const
      { return ((sx+tx-1)/tx) * ((sy+ty-1)/ty); }

    void chunk_info (int n, int &x0, int &x1, int &y0, int &y1) const
      {
      int ix = n%((sx+tx-1)/tx);
      int iy = n/((sx+tx-1)/tx);
      x0 = ix*tx; x1 = min(x0+tx, sx);
      y0 = iy*ty; y1 = min(y0+ty, sy);
      }
  };

void render (const vector<particle_sim> &p, arr2<COLOUR> &pic, 
      bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
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

#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      mpiMgr.allreduce_sum_raw
        (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres);
      if (mpiMgr.master())
        {
        if (a_eq_e)
          for(int ix=0;ix<xres;ix++)
            for(int iy=0;iy<yres;iy++)
              {
              pic[ix][iy].r=1-xexp(pic[ix][iy].r);
              pic[ix][iy].g=1-xexp(pic[ix][iy].g);
              pic[ix][iy].b=1-xexp(pic[ix][iy].b);
              }
        }
}

/////////////////////////////CUDA CODE///////////////////////////////////
#ifdef CUDA
void render_as_thread1 (const vector<particle_sim> &p, arr2<COLOUR> &pic, 
      bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
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

/*
It is no longer used since Dec 09. New one is render_as_thead1().
particle_splotch is not used any more.

void render_as_thread (const vector<particle_splotch> &p, arr2<COLOUR> &pic, 
      bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
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

		  COLOUR8 a=p[m].a, e, q;
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
*/

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
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
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

void render_cu_test1 (cu_particle_splotch *p, int n, cu_color **pic, 
      bool a_eq_e,double grayabsorb, void* buf)
      {
		cu_fragment_AeqE        *fbuf;
		cu_fragment_AneqE       *fbuf1;
		if (a_eq_e)
			fbuf =(cu_fragment_AeqE*) buf;
		else
			fbuf1 =(cu_fragment_AneqE*)buf;


      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
	  exptable xexp(-20.);

      int xres=800,  yres=800;
//      pic.fill(COLOUR(0,0,0));

#ifdef VS
	  work_distributor wd (xres,yres,xres,yres);
#else
      work_distributor wd (xres,yres,200,200);
#endif //ifdef VS

    #pragma omp parallel
{
      int chunk;
    #pragma omp for schedule(dynamic,1)
        int	fpos=0;
      for (chunk=0; chunk<wd.nchunks(); ++chunk)
        {
        int x0, x1, y0, y1;
        wd.chunk_info(chunk,x0,x1,y0,y1);
        arr2<COLOUR> lpic(x1-x0,y1-y0);
        lpic.fill(COLOUR(0,0,0));
        int x0s=x0, y0s=y0;
        x1-=x0; x0=0; y1-=y0; y0=0;


        for (unsigned int m=0; m<n; ++m)
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
                    fbuf[fpos].deltaR = (fac*a.r);
                    fbuf[fpos].deltaG = (fac*a.g);
                    fbuf[fpos].deltaB = (fac*a.b);
                  }
                else
                  {
//					  float tmp=xexp(fac*a.r);
//					  float	tmp1 =(1-q.r)*tmp;...
                  lpic[x][y].r = q.r+(lpic[x][y].r-q.r)*xexp(fac*a.r);
                  lpic[x][y].g = q.g+(lpic[x][y].g-q.g)*xexp(fac*a.g);
                  lpic[x][y].b = q.b+(lpic[x][y].b-q.b)*xexp(fac*a.b);
                  }//if a_eq_e
                }// if dsq<radsq
			  else
			  {
                if (a_eq_e)
                {
                    fbuf[fpos].deltaR =0.0;
                    fbuf[fpos].deltaG =0.0;
                    fbuf[fpos].deltaB =0.0;
                }
			  }
			  fpos++;
              }//y
            }//x
          }//for particle[m] 
        }//for this chunk
}//#pragma omp parallel

}

void render_cu_test (cu_particle_splotch *p, unsigned int size, arr2<COLOUR> &pic, 
      bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);

	  exptable xexp(MAX_EXP);
#ifdef CUDA_TEST_FRAGMENT
	  memset(fragBufWrittenByHost,0,sizeof(cu_fragment_AeqE)*30*100000*12);
#endif
      int xres = pic.size1(), yres=pic.size2();
      pic.fill(COLOUR(0,0,0));

	  work_distributor wd (xres,yres,xres,yres);

//estimate combination time
//float	t =0;
//counting fragments
long	fragments=0;
//long	PValid=0;

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


        for (unsigned int m=0; m<size; ++m)
          {
		  //as if it's on device
#ifdef CUDA_TEST_FRAGMENT
		  posFragBufH =p[m].posInFragBuf;
#endif
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
//to count valid particles with one chunk, debug only
//PValid++; 
//continue;
          COLOUR8 a(p[m].e.r, p[m].e.g, p[m].e.b) , e, q;
          if (!a_eq_e)
            {
            e.r=p[m].e.r; e.g=p[m].e.g; e.b=p[m].e.b;
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
#ifdef CUDA_TEST_EXP //for test exp table on device
	if (size_inout_buffer<100000)
	{
		inout_buffer[size_inout_buffer][0]=prefac1*dsq;
		inout_buffer[size_inout_buffer][1]=xexp(prefac1*dsq);
		size_inout_buffer++;
		continue;
	}
	else
		return;
#endif //if CUDA_TEST_EXP

                if (a_eq_e)
                  {
#ifdef CUDA_NOUSE
//estimate combination time, should remove later
//t =t+1.0;
//t =t+1.0;
//t =t+1.0; //3 times fload adding
//counting fragments
fragments++;
//continue;
#endif //CUDA_NOUSE

#ifdef CUDA_TEST_FRAGMENT
				  fragBufWrittenByHost[posFragBufH].deltaR =(fac*a.r);
				  fragBufWrittenByHost[posFragBufH].deltaG =(fac*a.g);
				  fragBufWrittenByHost[posFragBufH].deltaB =(fac*a.b);
#else
                  lpic[x][y].r += (fac*a.r);
                  lpic[x][y].g += (fac*a.g);
                  lpic[x][y].b += (fac*a.b);
#endif //CUDA_TEST_FRAGMENT
                  }
                else
                  {
                  lpic[x][y].r = q.r+(lpic[x][y].r-q.r)*xexp(fac*a.r);
                  lpic[x][y].g = q.g+(lpic[x][y].g-q.g)*xexp(fac*a.g);
                  lpic[x][y].b = q.b+(lpic[x][y].b-q.b)*xexp(fac*a.b);
                  }//if a_eq_e
                }// if dsq<radsq
#ifdef CUDA_TEST_FRAGMENT
				posFragBufH ++;				  	
#endif //CUDA_TEST_FRAGMENT
              }//y
            }//x
          }//for particle[m]
        for(int ix=0;ix<x1;ix++)
          for(int iy=0;iy<y1;iy++)
            pic[ix+x0s][iy+y0s]=lpic[ix][iy];
        }//for this chunk
//counting fragments, should remove later
//printf("\nfragments=%ld", fragments);
}//#pragma omp parallel

#ifdef CUDA_TEST_FRAGMENT
		return;				  	
#endif //CUDA_TEST_FRAGMENT

      mpiMgr.allreduce_sum_raw
        (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres);
      if (mpiMgr.master())
        {
        if (a_eq_e)
          for(int ix=0;ix<xres;ix++)
            for(int iy=0;iy<yres;iy++)
              {
              pic[ix][iy].r=1-xexp(pic[ix][iy].r);
              pic[ix][iy].g=1-xexp(pic[ix][iy].g);
              pic[ix][iy].b=1-xexp(pic[ix][iy].b);
              }
        }
}
#endif //ifdef CUDA

void add_colorbar(paramfile &params, arr2<COLOUR> &pic, vector<COLOURMAP> &amap)
{
  int xres = pic.size1(), yres=pic.size2();
  int offset=0;
  int ptypes = params.find<int>("ptypes",1);

  for(int itype=0;itype<ptypes;itype++)
    {
      if(params.find<bool>("color_is_vector"+dataToString(itype),false))
        {
           cout << " adding no color bar for type " << itype 
                << " as it is color vector ..." << endl;
        }
      else
        {
          cout << " adding color bar for type " << itype << " ..." << endl;
          for (int x=0; x<xres; x++)
            {  
               float64 temp=x/float64(xres);
               COLOUR e=amap[itype].Get_Colour(temp);
               for (int y=0; y<10; y++)
                 {
                   pic[x][yres-offset-1-y].r = e.r;
                   pic[x][yres-offset-1-y].g = e.g;
                   pic[x][yres-offset-1-y].b = e.b;
                 }
            }
          offset += 10;
        }
    }
}


void particle_normalize(paramfile &params, vector<particle_sim> &p, bool verbose) ///verbose means saying sth or not
{
  int ptypes = params.find<int>("ptypes",1);
  vector<bool> col_vector,log_int,log_col,asinh_col;
  vector<float32> mincol,maxcol,minint,maxint;
  col_vector.resize(ptypes);
  log_int.resize(ptypes);
  log_col.resize(ptypes);
  asinh_col.resize(ptypes);
  mincol.resize(ptypes);
  maxcol.resize(ptypes);
  minint.resize(ptypes);
  maxint.resize(ptypes);

  for(int itype=0;itype<ptypes;itype++)
    {
      log_int[itype] = params.find<bool>("intensity_log"+dataToString(itype),true);
      log_col[itype] = params.find<bool>("color_log"+dataToString(itype),true);
      asinh_col[itype] = params.find<bool>("color_asinh"+dataToString(itype),false);
      col_vector[itype] = params.find<bool>("color_is_vector"+dataToString(itype),false);
      mincol[itype]=1e30;
      maxcol[itype]=-1e30;
      minint[itype]=1e30;
      maxint[itype]=-1e30;
    }

  int npart=p.size();

  for (int m=0; m<npart; ++m) ///do log calcultions if demanded
    {
      if (log_int[p[m].type])
	p[m].I = log(p[m].I);
      get_minmax(minint[p[m].type], maxint[p[m].type], p[m].I);

      if (log_col[p[m].type])
	p[m].C1 = log(p[m].C1);
      if(asinh_col[p[m].type])
	p[m].C1 = my_asinh(p[m].C1);
      get_minmax(mincol[p[m].type], maxcol[p[m].type], p[m].C1);
      if (col_vector[p[m].type])
	{
	  if (log_col[p[m].type])
	    {
	      p[m].C2 = log(p[m].C2);
	      p[m].C3 = log(p[m].C3);
	    }   
	  if (asinh_col[p[m].type])
	    {
	      p[m].C2 = my_asinh(p[m].C2);
	      p[m].C3 = my_asinh(p[m].C3);
	    }
	  get_minmax(mincol[p[m].type], maxcol[p[m].type], p[m].C2);
	  get_minmax(mincol[p[m].type], maxcol[p[m].type], p[m].C3);
	}
    }


  for(int itype=0;itype<ptypes;itype++)
    {
      mpiMgr.allreduce_min(minint[itype]);
      mpiMgr.allreduce_min(mincol[itype]);
      mpiMgr.allreduce_max(maxint[itype]);
      mpiMgr.allreduce_max(maxcol[itype]);

      float minval_int = params.find<float>("intensity_	min"+dataToString(itype),minint[itype]);
      float maxval_int = params.find<float>("intensity_max"+dataToString(itype),maxint[itype]);
      float minval_col = params.find<float>("color_min"+dataToString(itype),mincol[itype]);
      float maxval_col = params.find<float>("color_max"+dataToString(itype),maxcol[itype]);

	  if(verbose && mpiMgr.master())
        {
           cout << " For particles of type " << itype << " : " << endl;
           cout << " From data: " << endl;
           cout << " Color Range:     " << mincol[itype] << " (min) , " <<
                                           maxcol[itype] << " (max) " << endl;
           cout << " Intensity Range: " << minint[itype] << " (min) , " <<
                                           maxint[itype] << " (max) " << endl;
           cout << " Restricted to: " << endl;
           cout << " Color Range:     " << minval_col << " (min) , " <<
                                           maxval_col << " (max) " << endl;
           cout << " Intensity Range: " << minval_int << " (min) , " <<
                                           maxval_int << " (max) " << endl;
        }

      for(int m=0; m<npart; ++m)
        {
          if(p[m].type == itype)///clamp into (min,max)
            {
               my_normalize(minval_int,maxval_int,p[m].I);
               my_normalize(minval_col,maxval_col,p[m].C1);
               if (col_vector[p[m].type])
	         {
	            my_normalize(minval_col,maxval_col,p[m].C2);
	            my_normalize(minval_col,maxval_col,p[m].C3);
	         }
            }
        }
    }
}


void particle_project(paramfile &params, vector<particle_sim> &p, VECTOR campos, VECTOR lookat, VECTOR sky)
{
  int res = params.find<int>("resolution",200);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);
  float64 xfac,dist;
  int npart=p.size();

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
  ///hereby matrix is setup

  for (long m=0; m<npart; ++m) ///now do the calculation
    {
      VECTOR v(p[m].x,p[m].y,p[m].z);
      v=trans.TransPoint(v);
      p[m].x=v.x; p[m].y=v.y; p[m].z=v.z;
    }

  bool projection = params.find<bool>("projection",true);

  if(!projection)
    {
      dist= (campos-lookat).Length();
      xfac=1./(fovfct*dist);
      cout << " Field of fiew: " << 1./xfac*2. << endl;
    }

  bool minhsmlpixel = params.find<bool>("minhsmlpixel",false);
     
  for (long m=0; m<npart; ++m) ///calculate ro, r
    {
      if(!projection)
        {
          p[m].x = res*.5 * (p[m].x+fovfct*dist)*xfac;
          p[m].y = res*.5 * (p[m].y+fovfct*dist)*xfac;
        }
      else
        {
          xfac=1./(fovfct*p[m].z);
          p[m].x = res*.5 * (p[m].x+fovfct*p[m].z)*xfac;
          p[m].y = res*.5 * (p[m].y+fovfct*p[m].z)*xfac;
        }
      p[m].ro = p[m].r;
      p[m].r = p[m].r *res*.5*xfac;
      if(minhsmlpixel)
         if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
	   {
             p[m].r = 0.5;
             p[m].ro = p[m].r/(res*.5*xfac);
	   }
    }
};

#ifdef CUDA
extern "C" 
void getCuTransformParams(cu_param_transform &para_trans,
paramfile &params, double c[3], double l[3], double s[3])
{
//	cout<<endl<<"\nRetrieve parameters for device transformation\n"<<endl;

	VECTOR	campos(c[0], c[1], c[2]), 
			lookat(l[0], l[1], l[2]), 
			sky(s[0], s[1], s[2]);
	
	int res = params.find<int>("resolution",200);
	double fov = params.find<double>("fov",45); //in degrees
	double fovfct = tan(fov*0.5*degr2rad);
	float64 xfac=0.0, dist=0.0;

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


void particle_colorize(paramfile &params, vector<particle_sim> &p, 
//                       vector<particle_splotch> &p2, 
                       vector<COLOURMAP> &amap, vector<COLOURMAP> &emap)
{
  int res = params.find<int>("resolution",200);
  int ycut0 = params.find<int>("ycut0",0);
  int ycut1 = params.find<int>("ycut1",res);
  float zmaxval = params.find<float>("zmax",1.e23);
  float zminval = params.find<float>("zmin",0.0);
  int ptypes = params.find<int>("ptypes",1);
  vector<bool> col_vector;
  vector<float64> brightness,grayabsorb;

  col_vector.resize(ptypes);
  brightness.resize(ptypes);
  grayabsorb.resize(ptypes);

//  p2.resize(0);
  for(int itype=0;itype<ptypes;itype++)
    {
      brightness[itype] = params.find<double>("brightness"+dataToString(itype),1.);
      grayabsorb[itype] = params.find<float>("gray_absorption"+dataToString(itype),0.2);
      col_vector[itype] = params.find<bool>("color_is_vector"+dataToString(itype),false);
    }
  float64 rfac=1.5;
  int npart=p.size();
  
  for (int m=0; m<npart; ++m)
    {
      p[m].active = 0;
      if (p[m].z<=0) continue;
      if (p[m].z<=zminval) continue;
      if (p[m].z>=zmaxval) continue;
      
      float64 r=p[m].r;
      float64 posx=p[m].x, posy=p[m].y;
      
      float64 rfacr=rfac*r;
      
      int minx=int(posx-rfacr+1);
      if (minx>=res) continue;
      minx=max(minx,0);
      int maxx=int(posx+rfacr+1);
      if (maxx<=0) continue;
      maxx=min(maxx,res);
      if (minx>=maxx) continue;
      int miny=int(posy-rfacr+1);
      if (miny>=ycut1) continue;
      miny=max(miny,ycut0);
      int maxy=int(posy+rfacr+1);
      if (maxy<=ycut0) continue;
      maxy=min(maxy,ycut1);
      if (miny>=maxy) continue;

      float64 col1=p[m].C1,col2=p[m].C2,col3=p[m].C3;
      clamp (0.0000001,0.9999999,col1);
      if (col_vector[p[m].type])
	{
          clamp (0.0000001,0.9999999,col2);
          clamp (0.0000001,0.9999999,col3);
	}
      float64 intensity=p[m].I;
      clamp (0.0000001,0.9999999,intensity);
      intensity *= brightness[p[m].type];
      
      COLOUR e;
      if (col_vector[p[m].type])
	{
          e.r=col1*intensity;
          e.g=col2*intensity;
          e.b=col3*intensity;
	}
      else
	e=amap[p[m].type].Get_Colour(col1)*intensity;

#ifdef CUDA_TEST_COLORMAP
//for CUDA TEST ONLY
	inout_buffer[size_inout_buffer][0] =p[m].type;
	inout_buffer[size_inout_buffer][1] =col1;
	inout_buffer[size_inout_buffer][2] =amap[p[m].type].Get_Colour(col1).r;
	inout_buffer[size_inout_buffer][3] =amap[p[m].type].Get_Colour(col1).g;
	inout_buffer[size_inout_buffer][4] =amap[p[m].type].Get_Colour(col1).b;
	size_inout_buffer++;
//CUDA test over
#endif
      COLOUR a=e;
      p[m].active = 1;
      p[m].e = e;
//      p2.push_back(particle_splotch(p[m].x, p[m].y,p[m].r,p[m].ro,e));
    }
}

void particle_sort(vector<particle_sim> &p, int sort_type, bool verbose)
{
  switch(sort_type)
    {
    case 0: 
      if(verbose && mpiMgr.master())
         cout << " skipped sorting ..." << endl;
      break;
    case 1: 
      if(verbose && mpiMgr.master())
         cout << " sorting by z ..." << endl;
      sort(p.begin(), p.end(), zcmp());
      break;
    case 2: 
      if(verbose && mpiMgr.master()) 
         cout << " sorting by value ..." << endl;
      sort(p.begin(), p.end(), vcmp1());
      break;
    case 3: if(verbose && mpiMgr.master()) 
         cout << " reverse sorting by value ..." << endl;
      sort(p.begin(), p.end(), vcmp2());
      break;
    case 4: 
      if(verbose && mpiMgr.master())
      cout << " sorting by size ..." << endl;
      sort(p.begin(), p.end(), hcmp());
      break;
    default:
      planck_fail("unknown sorting chosen ...");
      break;
    }
}


#endif

#ifdef INTERPOLATE

// Higher order interpolation would be:
// Time between snapshots (cosmology !)
//    dt=(z2t(h1.redshift)-z2t(h0.redshift))*0.7
// Velocity factors:
//    v_unit1=v_unit/l_unit/sqrt(h1.time)*dt
//    v_unit0=v_unit/l_unit/sqrt(h0.time)*dt
// Delta_v (cosmology)
//    vda=2*(x1-x0)-(v0*v_unit0+v1*v_unit1)
// Delta_t (0..1)
//     t=FLOAT(k)/FLOAT(nint) == frac (!)
// Interpolated positions: 
//    x=x0+v0*v_unit0*t+0.5*(v1*v_unit1-v0*v_unit0+vda)*t^2
// Interpolated velocities:
//    v=v0+t*(v1-v0)

void particle_interpolate(paramfile &params,
                          vector<particle_sim> &p,vector<particle_sim> &p1,
                          vector<particle_sim> &p2, double frac, double time1, 
                          double time2) 
{
  int i1=0,i2=0;

  cout << " Time1/2 = " << time1 << "," << time2 << endl;

#ifdef HIGH_ORDER_INTERPOLATION
  double h = params.find<double>("hubble",0.7);
  double O = params.find<double>("omega",0.3);
  double L = params.find<double>("lambda",0.7);
  double mparsck = 3.0856780e+24;
  double l_unit = params.find<double>("l_unit",3.0856780e+21);
  double v_unit = params.find<double>("v_unit",100000.00);
  double t1 = log(sqrt(L/O*time1*time1*time1)+sqrt((L/O*time1*time1*time1)+1))/1.5/sqrt(L)/h/1e7*mparsck;
  double t2 = log(sqrt(L/O*time2*time2*time2)+sqrt((L/O*time2*time2*time2)+1))/1.5/sqrt(L)/h/1e7*mparsck;
  double dt = (t2 - t1) * h;
  double v_unit1=v_unit/l_unit/sqrt(time1)*dt;
  double v_unit2=v_unit/l_unit/sqrt(time2)*dt;
  double vda_x,vda_y,vda_z;
#endif

  p.resize(0);
  while(i1 < p1.size() && i2 < p2.size())
    {
      int mode=-1;
      if(p1[i1].id == p2[i2].id)
         mode=0;
      if(p1[i1].id < p2[i2].id)
         mode=1;
      if(p1[i1].id > p2[i2].id)
         mode=2;
//      cout << "  (" << i1 << "," << i2 << "): " << p1[i1].id << "," << 
//              p2[i2].id << " => " << mode << " (" << p.size() << ")" << endl;
      switch(mode)
        {
         case 0:
                if (p1[i1].type != p2[i2].type)
                   planck_fail("interpolate: can not interpolate between different types !");
#ifdef HIGH_ORDER_INTERPOLATION
                vda_x = 2 * (p2[i2].x - p1[i1].x) - (p1[i1].vx * v_unit1 + p2[i2].vx * v_unit2);
                vda_y = 2 * (p2[i2].y - p1[i1].y) - (p1[i1].vy * v_unit1 + p2[i2].vy * v_unit2);
                vda_z = 2 * (p2[i2].z - p1[i1].z) - (p1[i1].vz * v_unit1 + p2[i2].vz * v_unit2);
#endif
                p.push_back(particle_sim(
#ifdef HIGH_ORDER_INTERPOLATION
                                         p1[i1].x + p1[i1].vx * v_unit1 * frac 
                                                  + 0.5 * (p2[i2].vx * v_unit2 - p1[i1].vx * v_unit1 + vda_x) * frac * frac,
                                         p1[i1].y + p1[i1].vy * v_unit1 * frac 
                                                  + 0.5 * (p2[i2].vy * v_unit2 - p1[i1].vy * v_unit1 + vda_y) * frac * frac,
                                         p1[i1].z + p1[i1].vz * v_unit1 * frac 
                                                  + 0.5 * (p2[i2].vz * v_unit2 - p1[i1].vz * v_unit1 + vda_z) * frac * frac,
#else
                                         (1-frac) * p1[i1].x  + frac*p2[i2].x,
                                         (1-frac) * p1[i1].y  + frac*p2[i2].y,
                                         (1-frac) * p1[i1].z  + frac*p2[i2].z,
#endif
                                         (1-frac) * p1[i1].r  + frac*p2[i2].r,
                                         (1-frac) * p1[i1].ro + frac*p2[i2].ro,
                                         (1-frac) * p1[i1].I  + frac*p2[i2].I,
                                         (1-frac) * p1[i1].C1 + frac*p2[i2].C1,
                                         (1-frac) * p1[i1].C2 + frac*p2[i2].C2,
                                         (1-frac) * p1[i1].C3 + frac*p2[i2].C3,
                                         p1[i1].type,p1[i1].active,p1[i1].e,p1[i1].id
#ifdef HIGH_ORDER_INTERPOLATION
                                        ,(1-frac) * p1[i1].vx  + frac*p2[i2].vx,
                                         (1-frac) * p1[i1].vy  + frac*p2[i2].vy,
                                         (1-frac) * p1[i1].vz  + frac*p2[i2].vz
#endif
                                         ));
                i1++;
                i2++;
                break;
         case 1:
                i1++;
                break;
         case 2:
                i2++;
                break;
         default:
                planck_fail("interpolate: should not happen !");
                break;
      }
    }
  if(mpiMgr.master())
    cout << " found " << p.size() << " common particles ..." << endl;

}
#endif
