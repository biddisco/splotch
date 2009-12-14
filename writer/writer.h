#ifndef WRITE_TGA_H
#define WRITE_TGA_H

#include<algorithm>
#ifdef VS
#include "cxxsupport/arr.h"
#include "cxxsupport/cxxutils.h"
#include "cxxsupport/paramfile.h"
#endif

#include "kernel/bstream.h"
#include "kernel/colour.h"

using namespace std;
using namespace RAYPP;


void write_tga(paramfile params, arr2<COLOUR> &pic, int res, string frame_name);

#endif
