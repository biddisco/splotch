#ifndef SPLOTCH_CUDA_H
#define SPLOTCH_CUDA_H

//data structs for using on device
//'d_' means device
struct d_particle_sim
{
	float x,y,z,r,ro,I,C1,C2,C3;
	int type;
#ifdef INTERPOLATE
	unsigned int id;
#endif
};

//functions
extern "C" void cu_init();
extern "C" void	cu_end();
extern "C" void	cu_range(d_particle_sim* p, unsigned int n);


#endif