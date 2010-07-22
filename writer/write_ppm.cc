#include <algorithm>

#include "kernel/bstream.h"
#include "writer/writer.h"

using namespace std;

void write_ppm(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
  {
  cout << " writing tga file '" << frame_name << "'" << endl;

  ofstream file(frame_name.c_str());

  file << "P3" << endl << pic.size1() << endl << pic.size2() << endl << 255 << endl; 
  for (tsize y=0; y<pic.size2(); ++y)
    {
      for (tsize x=0; x<pic.size1(); ++x)
	{
	  file << min(255,int(256*pic[x][y].r)) << " " 
	       << min(255,int(256*pic[x][y].g)) << " "
	       << min(255,int(256*pic[x][y].b)) << "   ";
	}
      file << endl;
    }
  }
