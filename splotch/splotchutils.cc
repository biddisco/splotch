#include "splotch/splotchutils.h"
#include "cxxsupport/mpi_support.h"
#include "kernel/transform.h"
#include "cuda/splotch_cuda.h"
using namespace std;

void render (const vector<particle_sim> &p, arr2<COLOUR> &pic, bool a_eq_e,
  double grayabsorb)
  {
  const float64 rfac=1.5;
  const float64 powtmp = pow(pi,1./3.);
  const float64 sigma0=powtmp/sqrt(2*pi);
  const float64 bfak=1./(2*sqrt(pi)*powtmp);
  exptable xexp(-20.);

  int xres = pic.size1(), yres=pic.size2();
  pic.fill(COLOUR(0,0,0));

  work_distributor wd (xres,yres,200,200);

#pragma omp parallel
{
  int chunk;
#pragma omp for schedule(dynamic,1)
  for (chunk=0; chunk<wd.nchunks(); ++chunk)
    {
    int x0, x1, y0, y1;
    wd.chunk_info(chunk,x0,x1,y0,y1);
    arr2<COLOUR8> lpic(x1-x0,y1-y0);
    arr<double> pre1(yres);
    lpic.fill(COLOUR8(0.,0.,0.));
    int x0s=x0, y0s=y0;
    x1-=x0; x0=0; y1-=y0; y0=0;

    for (unsigned int m=0; m<p.size(); ++m)
      if (p[m].active)
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

        COLOUR8 a=p[m].e, q(0,0,0);
        if (!a_eq_e)
          {
          COLOUR8 e=p[m].e;
          q=COLOUR8(e.r/(a.r+grayabsorb),e.g/(a.g+grayabsorb),e.b/(a.b+grayabsorb));
          }

        float64 radsq = rfacr*rfacr;
        float64 prefac1 = -0.5/(r*r*sigma0*sigma0);
        float64 prefac2 = -0.5*bfak/p[m].ro;
        for (int y=miny; y<maxy; ++y)
          pre1[y]=prefac2*xexp(prefac1*(y-posy)*(y-posy));

        for (int x=minx; x<maxx; ++x)
          {
          float64 xsq=(x-posx)*(x-posx);
          double pre2 = xexp(prefac1*xsq);

          for (int y=miny; y<maxy; ++y)
            {
            float64 dsq = (y-posy)*(y-posy) + xsq;
            if (dsq<radsq)
              {
              float64 fac = pre1[y]*pre2;
              if (a_eq_e)
                {
                lpic[x][y].r += (fac*a.r);
                lpic[x][y].g += (fac*a.g);
                lpic[x][y].b += (fac*a.b);
                }
              else
                {
                lpic[x][y].r += xexp.expm1(fac*a.r)*(lpic[x][y].r-q.r);
                lpic[x][y].g += xexp.expm1(fac*a.g)*(lpic[x][y].g-q.g);
                lpic[x][y].b += xexp.expm1(fac*a.b)*(lpic[x][y].b-q.b);
                } // if a_eq_e
              } // if dsq<radsq
            } // y
          } // x
        } // for particle[m]

    for (int ix=0;ix<x1;ix++)
      for (int iy=0;iy<y1;iy++)
        pic[ix+x0s][iy+y0s]=lpic[ix][iy];
    } // for this chunk
} // #pragma omp parallel

  mpiMgr.allreduceRaw
    (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres,MPI_Manager::Sum);
  if (mpiMgr.master() && a_eq_e)
    for (int ix=0;ix<xres;ix++)
      for (int iy=0;iy<yres;iy++)
        {
        pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
        pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
        pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
        }
  }

double my_asinh (double val)
  { return log(val+sqrt(1.+val*val)); }


void add_colorbar(paramfile &params, arr2<COLOUR> &pic, vector<COLOURMAP> &amap)
  {
  int xres = pic.size1(), yres=pic.size2();
  int offset=0;
  int ptypes = params.find<int>("ptypes",1);

  for(int itype=0;itype<ptypes;itype++)
    {
    if (params.find<bool>("color_is_vector"+dataToString(itype),false))
      cout << " adding no color bar for type " << itype
           << " as it is color vector ..." << endl;
    else
      {
      cout << " adding color bar for type " << itype << " ..." << endl;
      for (int x=0; x<xres; x++)
        {
        COLOUR e=amap[itype].getVal(x/float64(xres));
        for (int y=0; y<10; y++)
          pic[x][yres-offset-1-y] = e;
        }
      offset += 10;
      }
    }
  }


void particle_normalize(paramfile &params, vector<particle_sim> &p, bool verbose)
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

  for (int m=0; m<npart; ++m) //do log calculations if demanded
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

  for (long m=0; m<npart; ++m)
    {
    vec3 v(p[m].x,p[m].y,p[m].z);
    v=trans.TransPoint(v);
    p[m].x=v.x; p[m].y=v.y; p[m].z=v.z;
    }

  bool projection = params.find<bool>("projection",true);

  float64 dist = (campos-lookat).Length();
  float64 xfac = 1./(fovfct*dist);
  if (!projection)
    cout << " Field of fiew: " << 1./xfac*2. << endl;

  bool minhsmlpixel = params.find<bool>("minhsmlpixel",false);

  for (long m=0; m<npart; ++m)
    {
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

void particle_colorize(paramfile &params, vector<particle_sim> &p,
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
      e=amap[p[m].type].getVal(col1)*intensity;

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
    p[m].active = 1;
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

void particle_interpolate(paramfile &params, vector<particle_sim> &p,
  const vector<particle_sim> &p1, const vector<particle_sim> &p2, double frac,
  double time1, double time2)
  {
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
  tsize i1=0,i2=0;
  while(i1<p1.size() && i2<p2.size())
    {
    if (p1[i1].id==p2[i2].id)
      {
      planck_assert (p1[i1].type==p2[i2].type,
        "interpolate: can not interpolate between different types !");
#ifdef HIGH_ORDER_INTERPOLATION
      vda_x = 2 * (p2[i2].x-p1[i1].x) - (p1[i1].vx*v_unit1 + p2[i2].vx*v_unit2);
      vda_y = 2 * (p2[i2].y-p1[i1].y) - (p1[i1].vy*v_unit1 + p2[i2].vy*v_unit2);
      vda_z = 2 * (p2[i2].z-p1[i1].z) - (p1[i1].vz*v_unit1 + p2[i2].vz*v_unit2);
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
      }
    else if (p1[i1].id<p2[i2].id)
      i1++;
    else if (p1[i1].id>p2[i2].id)
      i2++;
    }

  if(mpiMgr.master())
    cout << " found " << p.size() << " common particles ..." << endl;
  }

#endif
