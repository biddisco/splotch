/*
 * Copyright (c) 2002-2005 Max-Planck-Society
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef PLANCK_CXXUTILS_H
#define PLANCK_CXXUTILS_H

#include <algorithm>
#include <string>
#include <map>
#include <cmath>
#include "message_error.h"
#include "constants.h"

//! Return \a true if | \a a-b | < \a epsilon * | \a b |, else \a false.
inline bool approx (double a, double b, double epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) < (epsilon*abs(b));
  }

//! Return \a true if | \a a-b | < \a epsilon, else \a false.
inline bool abs_approx (double a, double b, double epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) < epsilon;
  }

//! Returns the largest integer which is smaller than (or equal to) to \a arg.
inline int intfloor (double arg)
  {
  return (arg>=0) ? int(arg) : int(arg)-1;
  }

//! Returns the integer which is nearest to \a arg.
inline int planck_nint (double arg)
  {
  arg += 0.5;
  return (arg>=0) ? int(arg) : int(arg)-1;
  }

//! Returns the long integer which is nearest to \a arg.
inline long nlong (double arg)
  {
  arg += 0.5;
  return (arg>=0) ? long(arg) : long(arg)-1;
  }

//! Returns \a v1+v2 if \a v1<0, \a v1-v2 if \a v1>=v2, else \a v1.
/*! \a v1 can be positive or negative; \a v2 must be positive. */
template<typename T> inline T weak_modulo (T v1, T v2)
  { return (v1>=0) ? ((v1<v2) ? v1 : (v1-v2)) : (v1+v2); }

//! Returns the remainder of the division \a v1/v2.
/*! The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
inline double modulo (double v1, double v2)
  {
  using namespace std;
  return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
  }

//! Returns the remainder of the division \a v1/v2.
/*! The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
inline int modulo (int v1, int v2)
  { return (v1>=0) ? ((v1<v2) ? v1 : (v1%v2)) : ((v1%v2)+v2); }

//! Returns the remainder of the division \a v1/v2.
/*! The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
inline long modulo (long v1, long v2)
  { return (v1>=0) ? ((v1<v2) ? v1 : (v1%v2)) : ((v1%v2)+v2); }


//! Returns -1 if \a signvalue is negative, else +1
template<typename T> inline T sign (const T& signvalue)
  { return (signvalue>=0) ? 1 : -1; }

//! Returns the integer \a n, which fulfills \a n*n<=arg<(n+1)*(n+1).
inline unsigned int isqrt (unsigned int arg)
  {
  using namespace std;
  return unsigned (sqrt(arg+0.5));
  }

//! If the file \a filename is present, return \p true, else \p false.
bool file_present (const std::string &filename);

//! Checks the presence of the file \a filename
/*! If the file is not present, a Message_error is thrown. */
void assert_present (const std::string &filename);

//! Checks the absence of the file \a filename
/*! If the file is present, a Message_error is thrown. */
void assert_not_present (const std::string &filename);

//! Removes the file \a filename
void remove_file (const std::string &filename);

//! Test a condition, and throw a Message_error if it is false.
inline void planck_assert (bool testval, const std::string &msg)
  {
  if (testval) return;
  throw Message_error ("Assertion failed: "+msg);
  }
//! Test a condition, and throw a Message_error if it is false.
inline void planck_assert (bool testval, const char *msg)
  {
  if (testval) return;
  throw Message_error ("Assertion failed: "+std::string(msg));
  }


//! returns the string \a orig without leading and trailing whitespace
std::string trim (const std::string &orig);

//! returns a string containing the text representation of \a x.
template<typename T> std::string dataToString(const T &x);
template<> std::string dataToString (const bool &x);
template<> std::string dataToString (const std::string &x);
template<> std::string dataToString (const float &x);
template<> std::string dataToString (const double &x);

//! returns a string containing the text representation of \a x.
std::string intToString(int x, int width);

//! reads a value of a given datatype from a string
template<typename T> void stringToData (const std::string &x, T &value);
template<> void stringToData (const std::string &x, std::string &value);
template<> void stringToData (const std::string &x, bool &value);

//! reads a value of a given datatype from a string
template<typename T> inline T stringToData (const std::string &x)
  { T result; stringToData(x,result); return result; }

//! returns an index to the left of two interpolation values.
/*! \a begin points to an array containing a sequence of values
    sorted in ascending order. The length of the array is \a len.
    If \a val is lower than the first element, 0 is returned.
    If \a val is higher than the last element, \a len-2
    is returned. Else, the index of the largest element smaller
    than \a val is returned. */
template<typename T> inline int interpol_left
  (const T *begin, int len, const T &val)
  {
  const T *end = begin+len;
  const T *iter = std::lower_bound (begin, end, val);
  if (iter==begin) return 0;
  if (iter==end) return len-2;
  return (iter-begin)-1;
  }

//! returns an index to the nearest interpolation value.
/*! \a begin points to an array containing a sequence of values
    sorted in ascending order. The length of the array is \a len.
    If \a val is lower than the first element, 0 is returned.
    If \a val is higher than the last element, \a len-1 is returned.
    Else, the index of the nearest element within the sequence of
    values is returned. */
template<typename T> inline int interpol_nearest
  (const T *begin, int len, const T &val)
  {
  int left = interpol_left(begin, len, val);
  T delleft = val-(*(begin+left));
  T delright = (*(begin+left+1))-val;
  if (delright<0) return left+1;
  return (delright<delleft) ? (left+1) : left;
  }

//! returns \a atan2(y,x) if \a x!=0 or \a y!=0; else returns 0.
inline double safe_atan2 (double y, double x)
  {
  using namespace std;
  return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x);
  }

void announce_progress (int now, int total);
void announce_progress (double now, double last, double total);

//! Prints a banner containing \a name. Useful for displaying program names.
void announce (const std::string &name);


//! Parses the file \a filename and returns the key/value pairs in \a dict.
void parse_file (const std::string &filename,
  std::map<std::string,std::string> &dict);

//! returns an appropriate FITS repetition count for a map with \a npix pixels.
inline int healpix_repcount (int npix)
  {
  if (npix<1024) return 1;
  if ((npix%1024)==0) return 1024;
  return isqrt (npix/12);
  }
#endif
