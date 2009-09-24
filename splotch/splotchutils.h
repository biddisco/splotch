#ifndef SPLOTCHUTILS_H
#define SPLOTCHUTILS_H

#include "mpi_support.h"

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
  int type;
  };


struct particle_splotch
  {
  float32 x,y,r,ro;
  COLOUR a,e;

  particle_splotch (float32 x_, float32 y_, float32 r_, float32 ro_, const COLOUR &a_,
             const COLOUR &e_)
    : x(x_), y(y_), r(r_), ro(ro_), a(a_), e(e_) {}
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

class splotch_renderer
  {
  public:
    void render (const vector<particle_splotch> &p, arr2<COLOUR> &pic, 
                 bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
      exptable xexp(-20.);

      int xres = pic.size1(), yres=pic.size2();

      work_distributor wd (xres,yres,200,200);

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
                  }
                }
              }
            }
          }
        for(int ix=0;ix<x1;ix++)
          for(int iy=0;iy<y1;iy++)
            pic[ix+x0s][iy+y0s]=lpic[ix][iy];

        }
}

      mpiMgr.allreduce_sum(reinterpret_cast<float *>(&pic[0][0]),3*xres*yres);
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
  };

#endif
