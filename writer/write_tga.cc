#include <algorithm>

#include "kernel/bstream.h"
#include "writer/writer.h"

using namespace std;

void write_tga(paramfile &params, const arr2<COLOUR> &pic, tsize res,
  const string &frame_name)
  {
  tsize ycut0 = params.find<int>("ycut0",0);
  tsize ycut1 = params.find<int>("ycut1",res);

  cout << " writing tga file" << endl;
  tsize yres=ycut1-ycut0;
  const uint8 header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    res%256, res/256, yres%256, yres/256, 24, 32 };

  bofstream file(frame_name.c_str(),file_is_natural);

  file.put(&header[0],18);
  for (tsize y=ycut0; y<ycut1; ++y)
    {
    for (tsize x=0; x<res; ++x)
      {
      uint8 pix[3];
      pix[0] = uint8(min(255,int(256*pic[x][y].b)));
      pix[1] = uint8(min(255,int(256*pic[x][y].g)));
      pix[2] = uint8(min(255,int(256*pic[x][y].r)));
      file.put(&pix[0],3);
      }
    }
  }
