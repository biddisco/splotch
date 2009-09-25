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

void render (const vector<particle_splotch> &p, arr2<COLOUR> &pic, 
      bool a_eq_e,double grayabsorb)
      {
      const float64 rfac=1.5;
      const float64 powtmp = pow(Pi,1./3.);
      const float64 sigma0=powtmp/sqrt(2*Pi);
      const float64 bfak=1./(2*sqrt(Pi)*powtmp);
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


void particle_normalize(paramfile &params, vector<particle_sim> &p, bool verbose)
{
  int ptypes = params.find<int>("ptypes",1);
  bool col_vector[ptypes];
  bool log_int[ptypes];
  bool log_col[ptypes];
  bool asinh_col[ptypes];
  float32 mincol[ptypes], maxcol[ptypes], minint[ptypes], maxint[ptypes];

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

  for (int m=0; m<npart; ++m)
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

      float minval_int = params.find<float>("intensity_min"+dataToString(itype),minint[itype]);
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
          if(p[m].type == itype)
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


void paticle_project(paramfile &params, vector<particle_sim> &p, VECTOR campos, VECTOR lookat, VECTOR sky)
{
  int res = params.find<int>("resolution",200);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);
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
  
  for (long m=0; m<npart; ++m)
    {
      VECTOR v(p[m].x,p[m].y,p[m].z);
      v=trans.TransPoint(v);
      p[m].x=v.x; p[m].y=v.y; p[m].z=v.z;
    }
#ifdef PROJECTION_OFF
  float64 dist= (campos-lookat).Length();
  float64 xfac=1./(fovfct*dist);
  cout << " Field of fiew: " << 1./xfac*2. << endl;
#endif
      
  for (long m=0; m<npart; ++m)
    {
#ifdef PROJECTION_OFF
      p[m].x = res*.5 * (p[m].x+fovfct*dist)*xfac;
      p[m].y = res*.5 * (p[m].y+fovfct*dist)*xfac;
#else
      float64 xfac=1./(fovfct*p[m].z);
      p[m].x = res*.5 * (p[m].x+fovfct*p[m].z)*xfac;
      p[m].y = res*.5 * (p[m].y+fovfct*p[m].z)*xfac;
#endif
      p[m].ro = p[m].r;
      p[m].r = p[m].r *res*.5*xfac;
#ifdef MINHSMLPIXEL
      if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
	{
          p[m].r = 0.5;
          p[m].ro = p[m].r/(res*.5*xfac);
	}
#endif
    }
}


void particle_colorize(paramfile &params, vector<particle_sim> &p, 
                       vector<particle_splotch> &p2, 
                       vector<COLOURMAP> &amap, vector<COLOURMAP> &emap)
{
  int res = params.find<int>("resolution",200);
  int ycut0 = params.find<int>("ycut0",0);
  int ycut1 = params.find<int>("ycut1",res);
  float zmaxval = params.find<float>("zmax",1.e23);
  float zminval = params.find<float>("zmin",0.0);
  int ptypes = params.find<int>("ptypes",1);
  bool col_vector[ptypes];
  float64 brightness[ptypes];
  float64 grayabsorb[ptypes];
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
      
      COLOUR a=e;

      p2.push_back(particle_splotch(p[m].x, p[m].y,p[m].r,p[m].ro,a,e));
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
      planck_fail("unknown sorting choosen ...");
      break;
    }
}

#endif


