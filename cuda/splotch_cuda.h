#ifndef SPLOTCH_CUDA_H
#define SPLOTCH_CUDA_H

#include "cxxsupport\paramfile.h"

//data structs for using on device
//'d_' means device
struct cu_particle_sim
{
	float x,y,z,r,ro,I,C1,C2,C3;
	int type;
#ifdef INTERPOLATE
	unsigned int id;
#endif
};

#define MAX_P_TYPE 8 //('XXXX','TEMP','U','RHO','MACH','DTEG','DISS','VEL')
struct cu_param_range //parameters for range calculation
{
	int		ptypes;	//meaning how many types
	bool	col_vector[MAX_P_TYPE], log_int[MAX_P_TYPE], 
			log_col[MAX_P_TYPE],	asinh_col[MAX_P_TYPE];
	float	mincol[MAX_P_TYPE],		maxcol[MAX_P_TYPE],
			minint[MAX_P_TYPE],		maxint[MAX_P_TYPE];
};



//functions
extern "C" void cu_init();
extern "C" void	cu_end();
extern "C" void	cu_range(paramfile &params, cu_particle_sim* p, unsigned int n);


#endif