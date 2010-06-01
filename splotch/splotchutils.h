/*
 * Copyright (c) 2004-2010
 *              Martin Reinecke (1), Klaus Dolag (1)
 *               (1) Max-Planck-Institute for Astrophysics
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
#ifndef SPLOTCHUTILS_H
#define SPLOTCHUTILS_H

#include "cxxsupport/mpi_support.h"
#include "cxxsupport/paramfile.h"
#include "kernel/colour.h"
#include "cxxsupport/vec3.h"
#include "kernel/colourmap.h"
#include "cxxsupport/walltimer.h"

struct particle_sim
  {
  float32 x,y,z,r,I,C1,C2,C3;
  unsigned short type;
  bool active;
  COLOUR e;

#ifdef INTERPOLATE
  unsigned int id;
#ifdef HIGH_ORDER_INTERPOLATION
  float32 vx,vy,vz;
#endif
#endif
  particle_sim (float32 x_, float32 y_, float32 z_, float32 r_,
                float32 I_, float32 C1_, float32 C2_, float32 C3_, int type_,
                int active_, const COLOUR &e_
#ifdef INTERPOLATE
                , unsigned int id_
#ifdef HIGH_ORDER_INTERPOLATION
                , float32 vx_, float32 vy_, float32 vz_
#endif
#endif
                ): x(x_), y(y_), z(z_), r(r_), I(I_), C1(C1_), C2(C2_), C3(C3_), type(type_), active(active_), e(e_)
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

struct Normalizer
  {
  float32 minv, maxv;

  Normalizer ()
    : minv(1e37), maxv(-1e37) {}

  void collect (float32 val)
    {
    using namespace std;
    minv=min(minv,val); maxv=max(maxv,val);
    }

  void collect (const Normalizer &other)
    {
    using namespace std;
    minv=min(minv,other.minv); maxv=max(maxv,other.maxv);
    }

  void normAndClamp (float32 &val) const
    {
    using namespace std;
    val = (max(minv,min(maxv,val))-minv)/(maxv-minv);
    }
  };

class exptable
  {
  private:
    float32 expfac, taylorlimit;
    arr<float32> tab1, tab2;
    enum {
      nbits=10,
      dim1=1<<nbits,
      mask1=dim1-1,
      dim2=(1<<nbits)<<nbits,
      mask2=dim2-1,
      mask3=~mask2
      };

  public:
    exptable (float32 maxexp)
      : expfac(dim2/maxexp), tab1(dim1), tab2(dim1)
      {
      using namespace std;
      for (int m=0; m<dim1; ++m)
        {
        tab1[m]=exp(m*dim1/expfac);
        tab2[m]=exp(m/expfac);
        }
      taylorlimit = sqrt(2.f*abs(maxexp)/dim2);
      }

    float32 taylorLimit() const { return taylorlimit; }

    float32 operator() (float32 arg) const
      {
      int iarg= int(arg*expfac);
      if (iarg&mask3)
        return (iarg<0) ? 1.f : 0.f;
      return tab1[iarg>>nbits]*tab2[iarg&mask1];
      }
    float32 expm1(float32 arg) const
      {
      if (std::abs(arg)<taylorlimit) return arg;
      return operator()(arg)-1.f;
      }
  };

class work_distributor
  {
  private:
    int sx, sy, tx, ty, nx;

  public:
    work_distributor (int sx_, int sy_, int tx_, int ty_)
      : sx(sx_), sy(sy_), tx(tx_), ty(ty_), nx((sx+tx-1)/tx) {}

    int nchunks() const
      { return nx * ((sy+ty-1)/ty); }

    void chunk_info (int n, int &x0, int &x1, int &y0, int &y1) const
      {
      int ix = n%nx;
      int iy = n/nx;
      x0 = ix*tx; x1 = std::min(x0+tx, sx);
      y0 = iy*ty; y1 = std::min(y0+ty, sy);
      }
    void chunk_info_idx (int n, int &ix, int &iy) const
      {
      ix = n%nx;
      iy = n/nx;
      }
  };


void add_colorbar(paramfile &params, arr2<COLOUR> &pic,
  std::vector<COLOURMAP> &amap);

void timeReport();
void hostTimeReport(wallTimerSet &Timers);

void get_colourmaps (paramfile &params, std::vector<COLOURMAP> &amap);

double my_asinh (double val);

#endif // SPLOTCHUTILS_H
