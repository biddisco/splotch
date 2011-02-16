#include <algorithm>

#include "kernel/bstream.h"
#include "writer/writer.h"

using namespace std;

void write_bmp(paramfile &params, const arr2<COLOUR> &pic,
  const string &frame_name)
{
  cout << " writing bmp file '" << frame_name << "'" << endl;
  tsize xres=pic.size1(), yres=pic.size2();

  FILE *f;
  unsigned char *img = NULL;
  int filesize = 54 + 3*xres*yres;  
  img = (unsigned char *)malloc(3*xres*yres);
  memset(img,0,sizeof(img));
  int x,y,r,g,b;

  for(int i=0; i<xres; i++)
  {
      for(int j=0; j<yres; j++)
  {
      x=i; y=(yres-1)-j;
      r = pic[i][j].r*255;
      g = pic[i][j].g*255;
      b = pic[i][j].b*255;
      if (r > 255) r=255;
      if (g > 255) g=255;
      if (b > 255) b=255;
      img[(x+y*xres)*3+2] = (unsigned char)(r);
      img[(x+y*xres)*3+1] = (unsigned char)(g);
      img[(x+y*xres)*3+0] = (unsigned char)(b);
  }
  }

  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
  unsigned char bmppad[3] = {0,0,0};

  bmpfileheader[ 2] = (unsigned char)(filesize    );
  bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
  bmpfileheader[ 4] = (unsigned char)(filesize>>16);
  bmpfileheader[ 5] = (unsigned char)(filesize>>24);

  bmpinfoheader[ 4] = (unsigned char)(       xres    );
  bmpinfoheader[ 5] = (unsigned char)(       xres>> 8);
  bmpinfoheader[ 6] = (unsigned char)(       xres>>16);
  bmpinfoheader[ 7] = (unsigned char)(       xres>>24);
  bmpinfoheader[ 8] = (unsigned char)(       yres    );
  bmpinfoheader[ 9] = (unsigned char)(       yres>> 8);
  bmpinfoheader[10] = (unsigned char)(       yres>>16);
  bmpinfoheader[11] = (unsigned char)(       yres>>24);

  f = fopen(frame_name.c_str(),"wb");
  fwrite(bmpfileheader,1,14,f);
  fwrite(bmpinfoheader,1,40,f);
  for(int i=0; i<yres; i++)
  {
      fwrite(img+(xres*(yres-i-1)*3),3,xres,f);
      fwrite(bmppad,1,(4-(xres*3)%4)%4,f);
  }
  fclose(f);
}
