#ifndef SPLOTCHUTILS_CUDA_H
#define SPLOTCHUTILS_CUDA_H

#ifdef CUDA

#include <vector>
#include "cxxsupport/datatypes.h"
#include "kernel/colour.h"
#include "kernel/transform.h"
#include "splotch/splotchutils.h"
#include "cuda/splotch_cuda.h"

	//just for now that its still needed for cuda codes 2 Dec 2009
	struct particle_splotch
	{
		float32 x,y,r,ro;
		COLOUR	a, e;

		particle_splotch (float32 x_, float32 y_, float32 r_, float32 ro_, const COLOUR &a_, const COLOUR &e_)
		: x(x_), y(y_), r(r_), ro(ro_), a(a_), e(e_) {}
		particle_splotch () {}
	};
	
	//#define HOST_THREAD_RENDER
	#ifdef HOST_THREAD_RENDER
	struct particle_splotch;
	struct	param_render_thread
	{
		cu_particle_splotch	*p;
		int start, end;
	    bool a_eq_e;
		double grayabsorb;
		cu_color pic[][800];
//		void	**pic;
	};
	#endif //ifdef HOST_THREAD_RENDER
#endif //if CUDA

/////////////////////////////CUDA CODE///////////////////////////////////
#ifdef CUDA
void render_as_thread1 (const vector<particle_sim> &p, arr2<COLOUR> &pic, 
			bool a_eq_e,double grayabsorb);

#ifdef HOST_THREAD_RENDER
//DWORD WINAPI render_thread (const vector<particle_splotch> &p, int start, int end, arr2<COLOUR> &pic, 
//      bool a_eq_e,double grayabsorb)
DWORD WINAPI render_thread (param_render_thread *param);
#endif //ifdef HOST_THREAD_RENDER

extern "C" void getCuTransformParams(cu_param_transform &para_trans,
			  paramfile &params, double c[3], double l[3], double s[3]);

#endif	//CUDA

#endif
