#include <algorithm>

#include "kernel/bstream.h"
#include "writer/writer.h"

using namespace std;

void write_tga(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
  {
  cout << " writing tga file '" << frame_name << "'" << endl;
  tsize xres=pic.size1(), yres=pic.size2();
  const uint8 header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    xres%256, xres/256, yres%256, yres/256, 24, 32 };

  bofstream file(frame_name.c_str(),file_is_natural);

  file.put(&header[0],18);
  for (tsize y=0; y<xres; ++y)
    {
    for (tsize x=0; x<xres; ++x)
      {
      uint8 pix[3];
      pix[0] = uint8(min(255,int(256*pic[x][y].b)));
      pix[1] = uint8(min(255,int(256*pic[x][y].g)));
      pix[2] = uint8(min(255,int(256*pic[x][y].r)));
      file.put(&pix[0],3);
      }
    }
  }
