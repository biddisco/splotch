#ifndef SPLOTCHUTILS_H
#define SPLOTCHUTILS_H

#include "cxxsupport/mpi_support.h"
#include "cxxsupport/paramfile.h"
#include "kernel/colour.h"
#include "cxxsupport/vec3.h"
#include "kernel/colourmap.h"

using namespace std;

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
  unsigned int type:4;
  bool active:1;
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
  bool operator()(const particle_sim &p1, const particle_sim &p2) const
    { return p1.z>p2.z; }
  };

struct vcmp1
  {
  bool operator()(const particle_sim &p1, const particle_sim &p2) const
    { return p1.C1>p2.C1; }
  };

struct vcmp2
  {
  bool operator()(const particle_sim &p1, const particle_sim &p2) const
    { return p1.C1<p2.C1; }
  };

struct hcmp
  {
  bool operator()(const particle_sim &p1, const particle_sim &p2) const
    { return p1.r>p2.r; }
  };

template<typename T> void get_minmax (T &minv, T &maxv, T val)
  {
  minv=min(minv,val);
  maxv=max(maxv,val);
  }

template<typename T> void my_normalize (T minv, T maxv, T &val)
  { if (minv!=maxv) val=(val-minv)/(maxv-minv); }

template<typename T> void clamp (T minv, T maxv, T &val)
  { val = min(maxv, max(minv, val)); }


class exptable
  {
  private:
    float64 expfac, taylorlimit;
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
      using namespace std;
      for (int m=0; m<dim1; ++m)
        {
        tab1[m]=exp(m*dim1/expfac);
        tab2[m]=exp(m/expfac);
        }
      taylorlimit = sqrt(2.*abs(maxexp)/dim2);
      }

    float64 operator() (float64 arg) const
      {
      int iarg= int(arg*expfac);
      if (iarg&mask3)
        return (iarg<0) ? 1. : 0.;
      return tab1[iarg>>nbits]*tab2[iarg&mask1];
      }
    float64 expm1(float64 arg) const
      {
      if (abs(arg)<taylorlimit) return arg;
      return operator()(arg)-1.;
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
void particle_project(paramfile &params, vector<particle_sim> &p,
  const vec3 &campos, const vec3 &lookat, vec3 sky);
void particle_colorize(paramfile &params, vector<particle_sim> &p,
  vector<COLOURMAP> &amap, vector<COLOURMAP> &emap);
void particle_sort(vector<particle_sim> &p, int sort_type, bool verbose);
#ifdef INTERPOLATE
void particle_interpolate(paramfile &params,
                          vector<particle_sim> &p,vector<particle_sim> &p1,
                          vector<particle_sim> &p2, double frac, double time1,
                          double time2);
#endif

double my_asinh (double val);

#endif // SPLOTCHUTILS_H
