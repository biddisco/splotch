#include <fstream>
#include "splotch/splotchutils.h"
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/walltimer.h"

using namespace std;

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
    amap[itype].sortMap();
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
#ifndef CUDA
    cout << "Ranging Data (secs)        : " << wallTimers.acc("range") << endl;
    cout << "Transforming Data (secs)   : " << wallTimers.acc("transform") << endl;
    cout << "Sorting Data (secs)        : " << wallTimers.acc("sort") << endl;
    cout << "Coloring Sub-Data (secs)   : " << wallTimers.acc("coloring") << endl;
    cout << "Rendering Sub-Data (secs)  : " << wallTimers.acc("render") << endl;
#endif
    cout << "Postprocessing (secs)      : " << wallTimers.acc("postproc") << endl;
    cout << "Write Data (secs)          : " << wallTimers.acc("write") << endl;
    cout << "Total (secs)               : " << wallTimers.acc("full") << endl;
    wallTimers.start("full");
    }
  }
