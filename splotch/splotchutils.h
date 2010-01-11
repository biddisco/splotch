#ifndef SPLOTCHUTILS_H
#define SPLOTCHUTILS_H

#ifdef VS
#include "cxxsupport/mpi_support.h"
#endif
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/paramfile.h"
#include "kernel/colour.h"
#include "kernel/vector.h"
#include "kernel/colourmap.h"

#ifdef CUDA
#include "cuda/splotch_cuda.h"
#endif

using namespace std;

#ifdef CUDA
	//just for now that its still needed for cuda codes 2 Dec 2009
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
      int iarg= int(arg*expfac);
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
      bool a_eq_e,double grayabsorb);
void add_colorbar(paramfile &params, arr2<COLOUR> &pic, vector<COLOURMAP> &amap);
void particle_normalize(paramfile &params, vector<particle_sim> &p, bool verbose);
void particle_project(paramfile &params, vector<particle_sim> &p, VECTOR campos, VECTOR lookat, VECTOR sky);
void particle_colorize(paramfile &params, vector<particle_sim> &p, vector<COLOURMAP> &amap, vector<COLOURMAP> &emap);
void particle_sort(vector<particle_sim> &p, int sort_type, bool verbose);
#ifdef INTERPOLATE
void particle_interpolate(paramfile &params,
                          vector<particle_sim> &p,vector<particle_sim> &p1,
                          vector<particle_sim> &p2, double frac, double time1, 
                          double time2); 
#endif

double my_asinh (double val);

/////////////////////////////CUDA CODE///////////////////////////////////
#ifdef CUDA
void render_as_thread1 (const vector<particle_sim> &p, arr2<COLOUR> &pic, 
			bool a_eq_e,double grayabsorb);

#ifdef HOST_THREAD_RENDER
//DWORD WINAPI render_thread (const vector<particle_splotch> &p, int start, int end, arr2<COLOUR> &pic, 
//      bool a_eq_e,double grayabsorb)
DWORD WINAPI render_thread (param_render_thread *param);
#endif //ifdef HOST_THREAD_RENDER

extern "C" void getCuTransformParams(cu_param_transform &para_trans,
			  paramfile &params, double c[3], double l[3], double s[3]);

#endif	//CUDA


#endif            // SPLOTCHUTILS_H
