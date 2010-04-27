#ifndef SPLOTCH_HOST_H
#define SPLOTCH_HOST_H

#include "splotch/splotchutils.h"

void host_rendering(bool master, paramfile &params, 
                    std::vector<particle_sim> &particle_data, arr2<COLOUR> &pic,
                    vec3 &campos, vec3 &lookat, vec3 &sky, std::vector<COLOURMAP> &amap);

#endif
