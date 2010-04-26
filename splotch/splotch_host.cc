#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "splotch/splotch_host.h"
#include "kernel/transform.h"
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/walltimer.h"

using namespace std;

namespace {

void particle_normalize(paramfile &params, vector<particle_sim> &p, bool verbose)
  {
  int ptypes = params.find<int>("ptypes",1);
  arr<bool> col_vector(ptypes),log_int(ptypes),log_col(ptypes),asinh_col(ptypes);
  arr<float32> mincol(ptypes,1e30),maxcol(ptypes,-1e30),
               minint(ptypes,1e30),maxint(ptypes,-1e30);

  for(int itype=0;itype<ptypes;itype++)
    {
    log_int[itype] = params.find<bool>("intensity_log"+dataToString(itype),true);
    log_col[itype] = params.find<bool>("color_log"+dataToString(itype),true);
    asinh_col[itype] = params.find<bool>("color_asinh"+dataToString(itype),false);
    col_vector[itype] = params.find<bool>("color_is_vector"+dataToString(itype),false);
    }

  int npart=p.size();

  for (int m=0; m<npart; ++m) // do log calculations if requested
    {
    if (log_int[p[m].type])
      p[m].I = log10(p[m].I);
    get_minmax(minint[p[m].type], maxint[p[m].type], p[m].I);

    if (log_col[p[m].type])
      p[m].C1 = log10(p[m].C1);
    if (asinh_col[p[m].type])
      p[m].C1 = my_asinh(p[m].C1);
    get_minmax(mincol[p[m].type], maxcol[p[m].type], p[m].C1);
    if (col_vector[p[m].type])
      {
      if (log_col[p[m].type])
        {
        p[m].C2 = log10(p[m].C2);
        p[m].C3 = log10(p[m].C3);
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
    mpiMgr.allreduce(minint[itype],MPI_Manager::Min);
    mpiMgr.allreduce(mincol[itype],MPI_Manager::Min);
    mpiMgr.allreduce(maxint[itype],MPI_Manager::Max);
    mpiMgr.allreduce(maxcol[itype],MPI_Manager::Max);

    float minval_int = params.find<float>("intensity_min"+dataToString(itype),minint[itype]);
    float maxval_int = params.find<float>("intensity_max"+dataToString(itype),maxint[itype]);
    float minval_col = params.find<float>("color_min"+dataToString(itype),mincol[itype]);
    float maxval_col = params.find<float>("color_max"+dataToString(itype),maxcol[itype]);

    if (verbose && mpiMgr.master())
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


void particle_project(paramfile &params, vector<particle_sim> &p,
  const vec3 &campos, const vec3 &lookat, vec3 sky)
  {
  int res = params.find<int>("resolution",200);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);
  int npart=p.size();

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

  float64 dist = (campos-lookat).Length();
  float64 xfac = 1./(fovfct*dist);
  if (!projection)
    cout << " Field of fiew: " << 1./xfac*2. << endl;

  bool minhsmlpixel = params.find<bool>("minhsmlpixel",false);

#pragma omp parallel
{
  long m;
#pragma omp for schedule(static)
  for (m=0; m<npart; ++m)
    {
    vec3 v(p[m].x,p[m].y,p[m].z);
    v=trans.TransPoint(v);
    p[m].x=v.x; p[m].y=v.y; p[m].z=v.z;

    if (!projection)
      {
      p[m].x = res*.5*(p[m].x+fovfct*dist)*xfac;
      p[m].y = res*.5*(p[m].y+fovfct*dist)*xfac;
      }
    else
      {
      xfac=1./(fovfct*p[m].z);
      p[m].x = res*.5*(p[m].x+fovfct*p[m].z)*xfac;
      p[m].y = res*.5*(p[m].y+fovfct*p[m].z)*xfac;
      }
    p[m].ro = p[m].r;
    p[m].r = p[m].r *res*.5*xfac;
    if (minhsmlpixel && (p[m].r>=0.0))
      {
      p[m].r = sqrt(p[m].r*p[m].r + .5*.5);
      p[m].ro = p[m].r/(res*.5*xfac);
      }
    }
}
  }

void particle_colorize(paramfile &params, vector<particle_sim> &p,
  vector<COLOURMAP> &amap)
  {
  int res = params.find<int>("resolution",200);
  int ycut0 = params.find<int>("ycut0",0);
  int ycut1 = params.find<int>("ycut1",res);
  float zmaxval = params.find<float>("zmax",1.e23);
  float zminval = params.find<float>("zmin",0.0);
  int ptypes = params.find<int>("ptypes",1);
  arr<bool> col_vector(ptypes);
  arr<float64> brightness(ptypes),grayabsorb(ptypes);

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
    p[m].active = false;
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
      e=amap[p[m].type].getVal(col1)*intensity;

    p[m].active = true;
    p[m].e = e;
    }
  }

void particle_sort(vector<particle_sim> &p, int sort_type, bool verbose)
  {
  switch(sort_type)
    {
    case 0:
      if (verbose && mpiMgr.master())
        cout << " skipped sorting ..." << endl;
      break;
    case 1:
      if (verbose && mpiMgr.master())
        cout << " sorting by z ..." << endl;
      sort(p.begin(), p.end(), zcmp());
      break;
    case 2:
      if (verbose && mpiMgr.master())
        cout << " sorting by value ..." << endl;
      sort(p.begin(), p.end(), vcmp1());
      break;
    case 3:
      if (verbose && mpiMgr.master())
        cout << " reverse sorting by value ..." << endl;
      sort(p.begin(), p.end(), vcmp2());
      break;
    case 4:
      if (verbose && mpiMgr.master())
        cout << " sorting by size ..." << endl;
      sort(p.begin(), p.end(), hcmp());
      break;
    default:
      planck_fail("unknown sorting chosen ...");
      break;
    }
  }

const int chunkdim=100;

void render_new (vector<particle_sim> &p, arr2<COLOUR> &pic,
  bool a_eq_e, double grayabsorb, bool nopostproc)
  {
  planck_assert(a_eq_e || (mpiMgr.num_ranks()==1),
    "MPI only supported for A==E so far");

  int xres=pic.size1(), yres=pic.size2();
  int ncx=(xres+chunkdim-1)/chunkdim, ncy=(yres+chunkdim-1)/chunkdim;

  arr2<vector<uint32> > idx(ncx,ncy);
  double rcell=sqrt(2.)*(chunkdim*0.5-0.5);
  double cmid0=0.5*(chunkdim-1);

  const float64 rfac=1.5;
  const float64 powtmp = pow(pi,1./3.);
  const float64 sigma0=powtmp/sqrt(2*pi);
  const float64 bfak=1./(2*sqrt(pi)*powtmp);
  exptable xexp(-20.);

  pic.fill(COLOUR(0,0,0));

  for (tsize i=0; i<p.size(); ++i)
    {
    particle_sim &pp(p[i]);
    if (pp.active)
      {
      double pr = rfac*pp.r;
      pp.I = -0.5/(pp.r*pp.r*sigma0*sigma0);
      pp.r=pr;
      if (!a_eq_e)
        {
        pp.C1=pp.e.r/(pp.e.r+grayabsorb);
        pp.C2=pp.e.g/(pp.e.g+grayabsorb);
        pp.C3=pp.e.b/(pp.e.b+grayabsorb);
        }
      float64 intens = -0.5*bfak/pp.ro;
      pp.e.r*=intens; pp.e.g*=intens; pp.e.b*=intens;

      int minx=max(0,int(pp.x-pr+1)/chunkdim);
      int maxx=min(ncx-1,int(pp.x+pr)/chunkdim);
      int miny=max(0,int(pp.y-pr+1)/chunkdim);
      int maxy=min(ncy-1,int(pp.y+pr)/chunkdim);
      double px=pp.x, py=pp.y;
      double sumsq=(rcell+pr)*(rcell+pr);
      for (int ix=minx; ix<=maxx; ++ix)
        {
        double cx=cmid0+ix*chunkdim;
        for (int iy=miny; iy<=maxy; ++iy)
          {
          double cy=cmid0+iy*chunkdim;
          double rtot2 = (px-cx)*(px-cx) + (py-cy)*(py-cy);
          if (rtot2<sumsq)
            idx[ix][iy].push_back(i);
          }
        }
      }
    }

  work_distributor wd (xres,yres,chunkdim,chunkdim);
#pragma omp parallel
{
  arr<float64> pre1(chunkdim);
  arr2<COLOUR8> lpic(chunkdim,chunkdim);
  int chunk;
#pragma omp for schedule(dynamic,1)
  for (chunk=0; chunk<wd.nchunks(); ++chunk)
    {
    int x0, x1, y0, y1;
    wd.chunk_info(chunk,x0,x1,y0,y1);
    lpic.fast_alloc(x1-x0,y1-y0);
    lpic.fill(COLOUR8(0,0,0));
    int x0s=x0, y0s=y0;
    x1-=x0; x0=0; y1-=y0; y0=0;
    int cx, cy;
    wd.chunk_info_idx(chunk,cx,cy);
    const vector<uint32> &v(idx[cx][cy]);

    for (tsize m=0; m<v.size(); ++m)
      {
      const particle_sim &pp(p[v[m]]);
      float64 r=pp.r;
      float64 posx=pp.x, posy=pp.y;
      posx-=x0s; posy-=y0s;
      int minx=int(posx-r+1);
      minx=max(minx,x0);
      int maxx=int(posx+r+1);
      maxx=min(maxx,x1);
      int miny=int(posy-r+1);
      miny=max(miny,y0);
      int maxy=int(posy+r+1);
      maxy=min(maxy,y1);

      COLOUR8 a=pp.e;

      float64 radsq = r*r;
      float64 stp = pp.I;
      float64 att_max=xexp(stp*radsq);
      for (int y=miny; y<maxy; ++y)
        pre1[y]=xexp(stp*(y-posy)*(y-posy));

      if (a_eq_e)
        {
        for (int x=minx; x<maxx; ++x)
          {
          double pre2 = xexp(stp*(x-posx)*(x-posx));
          for (int y=miny; y<maxy; ++y)
            {
            float64 att = pre1[y]*pre2;
//FIXME: this if() seems to hurt performance
            if (att>att_max)
              {
              lpic[x][y].r += att*a.r;
              lpic[x][y].g += att*a.g;
              lpic[x][y].b += att*a.b;
              }
            }
          }
        }
      else
        {
        COLOUR8 q(pp.C1,pp.C2,pp.C3);

        for (int x=minx; x<maxx; ++x)
          {
          double pre2 = xexp(stp*(x-posx)*(x-posx));
          for (int y=miny; y<maxy; ++y)
            {
            float64 att = pre1[y]*pre2;
//FIXME: this if() seems to hurt performance
            if (att>att_max)
              {
              lpic[x][y].r += xexp.expm1(att*a.r)*(lpic[x][y].r-q.r);
              lpic[x][y].g += xexp.expm1(att*a.g)*(lpic[x][y].g-q.g);
              lpic[x][y].b += xexp.expm1(att*a.b)*(lpic[x][y].b-q.b);
              }
            }
          }
        }
      } // for particle

    for (int ix=0;ix<x1;ix++)
      for (int iy=0;iy<y1;iy++)
        pic[ix+x0s][iy+y0s]=lpic[ix][iy];
    } // for this chunk
} // #pragma omp parallel

  mpiMgr.allreduceRaw
    (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres,MPI_Manager::Sum);

  if (!nopostproc)
    if (mpiMgr.master() && a_eq_e)
      for (int ix=0;ix<xres;ix++)
        for (int iy=0;iy<yres;iy++)
          {
          pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
          pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
          pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
          }
  }


} // unnamed namespace

void host_processing(bool master, paramfile &params, long npart_all, 
                     vector<particle_sim> &particles,
                    vec3 &campos, vec3 &lookat, vec3 &sky, vector<COLOURMAP> &amap)
  {
// -----------------------------------
// ----------- Ranging ---------------
// -----------------------------------
  wallTimers.start("range");
  if (master)
    cout << endl << "host: ranging values (" << npart_all << ") ..." << endl;
  particle_normalize(params,particles,true); ///does log calculations and clamps data
  wallTimers.stop("range");

// -------------------------------------
// ----------- Transforming ------------
// -------------------------------------
  wallTimers.start("transform");
  if (master)
    cout << endl << "host: applying geometry (" << npart_all << ") ..." << endl;
  particle_project(params, particles, campos, lookat, sky);
  wallTimers.stop("transform");

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
  wallTimers.start("sort");
  if (master)
    (mpiMgr.num_ranks()>1) ?
      cout << endl << "host: applying local sort ..." << endl :
      cout << endl << "host: applying sort (" << particles.size() << ") ..." << endl;
  int sort_type = params.find<int>("sort_type",1);
  particle_sort(particles,sort_type,true);
  wallTimers.stop("sort");

// ------------------------------------
// ----------- Coloring ---------------
// ------------------------------------
  wallTimers.start("coloring");
  if (master)
    cout << endl << "host: calculating colors (" << npart_all << ") ..." << endl;
  particle_colorize(params, particles, amap);
  wallTimers.stop("coloring");

  }

