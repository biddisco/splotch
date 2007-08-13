/*
 *  Ray++ - Object-oriented ray tracing library
 *  Copyright (C) 1998-2001 Martin Reinecke and others.
 *  See the AUTHORS file for more information.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  See the README file for more information.
 */

namespace RAYPP {

class TYPESELECT
  {
  private:
    template<bool a, bool b, bool c, bool d,
             typename A, typename B, typename C, typename D>
      struct SELECT4
      {
      typedef void RET;
      };
    template<bool b, bool c, bool d,
             typename A, typename B, typename C, typename D>
      struct SELECT4<true, b, c, d, A, B, C, D>
      {
      typedef A RET;
      };
    template<bool c, bool d, typename A, typename B, typename C, typename D>
      struct SELECT4<false, true, c, d, A, B, C, D>
      {
      typedef B RET;
      };
    template<bool d, typename A, typename B, typename C, typename D>
      struct SELECT4<false, false, true, d, A, B, C, D>
      {
      typedef C RET;
      };
    template<typename A, typename B, typename C, typename D>
      struct SELECT4<false, false, false, true, A, B, C, D>
      {
      typedef D RET;
      };

  public:
    template<int num_bytes> struct GET_INT
      {
      typedef typename
        SELECT4 <(sizeof(signed char) >= num_bytes),
                 (sizeof(short int  ) >= num_bytes),
                 (sizeof(int        ) >= num_bytes),
                 (sizeof(long int   ) >= num_bytes),
                 signed char, short int, int, long int>::RET RET;
      };
    template<int num_bytes> struct GET_UINT
      {
      typedef typename
        SELECT4 <(sizeof(unsigned  char    ) >= num_bytes),
                 (sizeof(unsigned short int) >= num_bytes),
                 (sizeof(unsigned       int) >= num_bytes),
                 (sizeof(unsigned  long int) >= num_bytes),
                  unsigned char, unsigned short int,
                  unsigned int, unsigned long int>::RET RET;
      };

    template<int num_bytes> struct GET_FLOAT
      {
      typedef typename
        SELECT4 <(sizeof(float      ) >= num_bytes),
                 (sizeof(double     ) >= num_bytes),
                 (sizeof(long double) >= num_bytes),
                 true, float, double, long double, void>::RET RET;
      };
  };

} //namespace RAYPP
