#ifndef READER_H
#define READER_H

#include "splotch/splotchutils.h"

void gadget_reader(paramfile &params, vector<particle_sim> &p, int snr, double *time);
void gadget_millenium_reader(paramfile &params, vector<particle_sim> &p, int snr, double *time);
long bin_reader_tab (paramfile &params, vector<particle_sim> &points, float *maxr, float *minr,int mype, int npes);
long bin_reader_block (paramfile &params, vector<particle_sim> &points, float *maxr, float *minr, int mype, int npes);
long bin_reader_block_mpi (paramfile &params, vector<particle_sim> &points, float *maxr, float *minr, int mype, int npes);
#endif
