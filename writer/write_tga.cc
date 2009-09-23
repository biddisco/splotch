#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>
#include "arr.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"

using namespace std;
using namespace RAYPP;


void write_tga(paramfile params, arr2<COLOUR> &pic, int res, string frame_name)
{
  int ycut0 = params.find<int>("ycut0",0);
  int ycut1 = params.find<int>("ycut1",res);

  cout << "writing" << endl;
  int yres=ycut1-ycut0;
  const byte header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    res%256, res/256, yres%256, yres/256, 24, 32 };

  bofstream file(frame_name.c_str(),file_is_natural);

  for (int m=0; m<18; ++m) file << header[m];
  for (int y=ycut0; y<ycut1; ++y)
    {
    for (int x=0; x<res; ++x)
      {

      pic[x][y].b=min(float64(1.), max(float64(0.), float64(pic[x][y].b)));
      byte pix = byte(min(255,int(256*pic[x][y].b)));
      // patch // byte pix = min(byte(255),byte(256*pic[x][y].b));
      file << pix;
      pic[x][y].g=min(float64(1.), max(float64(0.), float64(pic[x][y].g)));
      pix = byte(min(255,int(256*pic[x][y].g)));
      // patch // pix = min(byte(255),byte(256*pic[x][y].g));
      file << pix;
      pic[x][y].r=min(float64(1.), max(float64(0.), float64(pic[x][y].r)));
      pix = byte(min(255,int(256*pic[x][y].r)));
      // patch // pix = min(byte(255),byte(256*pic[x][y].r));
      file << pix;
      }
    }

// END OF FUNCTION
  }
