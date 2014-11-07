/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Classes for creation and output of image files
 *
 *  Copyright (C) 2003-2012 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <fstream>
#include <sstream>
#include "ls_image.h"
#include "bstream.h"
#include "font_data.inc"
#include "string_utils.h"
#include "share_utils.h"
#include <stdio.h>
#include <string.h>

using namespace std;

const MP_Font medium_bold_font = { 0, 128, 7, 13, medium_bold_font_data };
const MP_Font giant_font = { 0, 128, 9, 15, giant_font_data };

void Palette::setPredefined (int num)
  {
  fv.clear(); cv.clear();
  switch(num)
    {
    case 0:
      add(0,Colour(0,0,0));
      add(1,Colour(1,1,1));
      break;
    case 1:
      add(0,Colour(0,0,0));
      add(0.4f,Colour(0,0,0.5f));
      add(0.75f,Colour(0,0.6f,1));
      add(1,Colour(1,1,1));
      break;
    case 4:
      add(0,Colour(0,0,.5f));
      add(0.15f,Colour(0,0,1));
      add(0.4f,Colour(0,1,1));
      add(0.7f,Colour(1,1,0));
      add(0.9f,Colour(1,.33f,0));
      add(1,Colour(.5f,0,0));
      break;
    default:
      planck_fail("Palette #"+dataToString(num)+" not yet supported.");
    }
  }

void LS_Image::write_char (int xpos, int ypos, const Colour &col, char c,
  int scale)
  {
  planck_assert ((c>=font.offset) && (c<font.offset+font.num_chars),
    "write_char: character out of range");
  for (int i=0; i<font.xpix; ++i)
    for (int j=0; j<font.ypix; ++j)
      {
      int ofs = (c-font.offset)*font.xpix*font.ypix + j*font.xpix + i;
      if (font.data[ofs]!=' ')
        for (int m=0; m<scale; ++m)
          for (int n=0; n<scale; ++n)
            put_pixel(xpos+scale*i+m,ypos+scale*j+n,col);
      }
  }

LS_Image::LS_Image ()
  : font(medium_bold_font) {}

LS_Image::LS_Image (int xres, int yres)
  : font(medium_bold_font), pixel(xres,yres,Colour(0,0,0)) {}

void LS_Image::annotate (int xpos, int ypos, const Colour &col,
  const string &text, int scale)
  {
  for (tsize m=0; m<text.length(); ++m)
    write_char(xpos+m*scale*font.xpix, ypos, col, text[m],scale);
  }

void LS_Image::annotate_centered (int xpos, int ypos, const Colour &col,
  const string &text, int scale)
  {
  xpos-=(scale*text.length()*font.xpix)/2;
  ypos-=scale*font.ypix/2;
  annotate (xpos,ypos,col,text,scale);
  }

void LS_Image::set_font (const MP_Font &fnt)
  { font = fnt; }

void LS_Image::write_TGA (const string &file) const
  {
  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file '" + file + "'");

  tsize xres=pixel.size1(), yres=pixel.size2();

  bostream bo(out);
  const uint8 header[18] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    uint8(xres%256), uint8(xres/256), uint8(yres%256), uint8(yres/256), 24, 32};

  bo.put (header, 18);

  for (tsize j=0; j<yres; ++j)
    for (tsize i=0; i<xres; ++i)
      bo << pixel[i][j].b << pixel[i][j].g << pixel[i][j].r;

  planck_assert(out,"error writing output file '" + file + "'");
  }

namespace {

void write_equal_range (const arr<Colour8> &px, tsize begin, tsize end,
  bostream &file)
  {
  chunkMaker cm (end-begin,128);
  uint64 cbeg, csz;
  while (cm.getNext(cbeg,csz))
    {
    file << uint8(csz-1+128);
    file << px[begin].b << px[begin].g << px[begin].r;
    }
  }
void write_unequal_range (const arr<Colour8> &px, tsize begin, tsize end,
  bostream &file)
  {
  chunkMaker cm (end-begin,128);
  uint64 cbeg, csz;
  while (cm.getNext(cbeg,csz))
    {
    file << uint8(csz-1);
    for (tsize cnt=begin+cbeg; cnt< begin+cbeg+csz; ++cnt)
      file << px[cnt].b << px[cnt].g << px[cnt].r;
    }
  }

} // unnamed namespace

void LS_Image::write_TGA_rle(const string &file) const
  {
  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file '" + file + "'");

  tsize xres=pixel.size1(), yres=pixel.size2();

  bostream bo(out);
  const uint8 header[18] = { 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    uint8(xres%256), uint8(xres/256), uint8(yres%256), uint8(yres/256), 24, 32};

  bo.put(header,18);
  for (tsize y=0; y<yres; ++y)
    {
    arr<Colour8> px(xres);
    for (tsize x=0; x<xres; ++x) px[x] = pixel[x][y];
    tsize xstart=0;
    while (xstart<xres)
      {
      if (xstart==xres-1)
        {
        write_unequal_range (px,xstart,xstart+1,bo);
        xstart=xres;
        }
      else
        {
        if (px[xstart+1]==px[xstart]) // range of equal pixels
          {
          tsize xend=xstart+2;
          while ((xend<xres) && (px[xend]==px[xstart])) ++xend;
          write_equal_range (px,xstart,xend,bo);
          xstart=xend;
          }
        else
          {
          tsize xend=xstart+2;
          while ((xend<xres) && (px[xend]!=px[xend-1])) ++xend;
          write_unequal_range (px,xstart,xend,bo);
          xstart=xend;
          }
        }
      }
    }
  planck_assert(out,"error writing output file '" + file + "'");
  }

void LS_Image::write_PPM (const string &file) const
  {
  ofstream out(file.c_str(), ios_base::out | ios_base::binary);
  planck_assert(out, "could not create file '" + file + "'");

  tsize xres=pixel.size1(), yres=pixel.size2();

  bostream bo(out);

  ostringstream header;
  header << "P6" << endl << xres << endl << yres << endl << 255 << endl;
  string hdrdata = header.str();
  bo.put(hdrdata.c_str(),hdrdata.size());

  for (tsize j=0; j<yres; ++j)
    for (tsize i=0; i<xres; ++i)
      bo << pixel[i][j].r << pixel[i][j].g << pixel[i][j].b;

  planck_assert(out,"error writing output file '" + file + "'");
  }

void LS_Image::write_bmp (const string &file) const
{
  cout << " writing bmp file '" << file.c_str() << "'" << endl;
  tsize xres=pixel.size1(), yres=pixel.size2();

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
      x=i; // y=(yres-1)-j;
      y = j;
      r = pixel[i][j].r;
      g = pixel[i][j].g;
      b = pixel[i][j].b;
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

  f = fopen(file.c_str(),"wb");
  fwrite(bmpfileheader,1,14,f);
  fwrite(bmpinfoheader,1,40,f);
  for(int i=0; i<yres; i++)
  {
      fwrite(img+(xres*(yres-i-1)*3),3,xres,f);
      fwrite(bmppad,1,(4-(xres*3)%4)%4,f);
  }
  fclose(f);
}
