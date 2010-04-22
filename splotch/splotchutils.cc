#include <fstream>
#include "splotch/splotchutils.h"
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/walltimer.h"

using namespace std;

void create_new_particles (const vector<particle_sim> &in, bool a_eq_e,
  float64 gray, vector<particle_new> &out, vector<locinfo> &loc, vector<COLOUR> &qvec)
  {
  const float64 powtmp = pow(pi,1./3.);
  const float64 sigma0=powtmp/sqrt(2*pi);
  const float64 bfak=1./(2*sqrt(pi)*powtmp);
  const int maxpix=65000;

  tsize nactive=0;
  for (tsize i=0; i<in.size(); ++i)
    if (in[i].active) ++nactive;

  out.reserve(nactive);
  loc.reserve(nactive);
  if (!a_eq_e) qvec.reserve(nactive);
  for (tsize i=0; i< in.size(); ++i)
    {
    if (in[i].active)
      {
      particle_new p;
      p.x = in[i].x;
      p.y = in[i].y;
      p.a = in[i].e;
      if (!a_eq_e)
        qvec.push_back (COLOUR (p.a.r/(p.a.r+gray),p.a.g/(p.a.g+gray),p.a.b/(p.a.b+gray)));

      p.a = p.a *(-0.5*bfak/in[i].ro);
      const float64 min_change=8e-5;
      p.steepness = -0.5/(in[i].r*in[i].r*sigma0*sigma0);
//      float64 amax=max(abs(p.a.r),max(abs(p.a.g),abs(p.a.b)));
//      const float64 min_change=8e-5;
//      float64 attenuation = min_change/amax;
      float64 attenuation = 3.7e-2; // equivalent to rfac=1.5
      p.rmax = sqrt(max(0.,log(attenuation)/p.steepness));
//p.rmax=3*in[i].r;
      out.push_back(p);

      locinfo l;
      l.minx = uint16(max(0,min(maxpix,int(p.x-p.rmax+1))));
      l.maxx = uint16(max(0,min(maxpix,int(p.x+p.rmax+1))));
      l.miny = uint16(max(0,min(maxpix,int(p.y-p.rmax+1))));
      l.maxy = uint16(max(0,min(maxpix,int(p.y+p.rmax+1))));
      loc.push_back(l);
      }
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


void get_colourmaps (paramfile &params, vector<COLOURMAP> &amap)
  {
  int ptypes = params.find<int>("ptypes",1);
  bool master = mpiMgr.master();
  amap.resize(ptypes);

  if (master)
    cout << "building color maps (" << ptypes << ")..." << endl;
  for (int itype=0;itype<ptypes;itype++)
    {
    if (params.find<bool>("color_is_vector"+dataToString(itype),false))
      {
      if (master)
        cout << " color of ptype " << itype << " is vector, so no colormap to load ..." << endl;
      }
    else
      {
      ifstream infile (params.find<string>("palette"+dataToString(itype)).c_str());
      planck_assert (infile,"could not open palette file  <" +
        params.find<string>("palette"+dataToString(itype)) + ">");
      string dummy;
      int nColours;
      infile >> dummy >> dummy >> nColours;
      if (master)
        cout << " loading " << nColours << " entries of color table of ptype " << itype << endl;
      double step = 1./(nColours-1);
      for (int i=0; i<nColours; i++)
        {
        float rrr,ggg,bbb;
        infile >> rrr >> ggg >> bbb;
        amap[itype].addVal(i*step,COLOUR(rrr/255,ggg/255,bbb/255));
        }
      }
    }
  }

void timeReport()
  {
  if (mpiMgr.master())
    {
    wallTimers.stop("full");
    cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "Setup Data (secs)          : " << wallTimers.acc("setup") << endl;
    cout << "Read Data (secs)           : " << wallTimers.acc("read") << endl;
    cout << "Ranging Data (secs)        : " << wallTimers.acc("range") << endl;
    cout << "Transforming Data (secs)   : " << wallTimers.acc("transform") << endl;
    cout << "Sorting Data (secs)        : " << wallTimers.acc("sort") << endl;
    cout << "Coloring Sub-Data (secs)   : " << wallTimers.acc("coloring") << endl;
    cout << "Rendering Sub-Data (secs)  : " << wallTimers.acc("render") << endl;
    cout << "Write Data (secs)          : " << wallTimers.acc("write") << endl;
    cout << "Total (secs)               : " << wallTimers.acc("full") << endl;
    wallTimers.start("full");
    }
  }
