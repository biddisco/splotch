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
	
	//#define CUDA_TEST_FRAGMENT
	
	#ifdef CUDA_TEST_FRAGMENT
		//the fragment buffer at host
		unsigned long	posFragBufH=0;
		cu_fragment_AeqE	fragBuf[30*100000*12];//for current test only! 10-29
		cu_fragment_AeqE	fragBufWrittenByHost[30*100000*12];//for current test only! 10-29
	#endif //ifdef CUDA_TEST_FRAGMENT

	//a array used for debug only
	#ifdef CUDA_TEST_COLORMAP
		int		size_inout_buffer=0;	//for debug only.
		float	inout_buffer[500000][5];//for debug only. ptype,input,r,g,b, 
	#endif //if CUDA_TEST_COLORMAP

	//#define CUDA_TEST_EXP
	#ifdef CUDA_TEST_EXP
		int		size_inout_buffer=0;	//for debug only.
		float	inout_buffer[100000][2];//for debug only. in, out
	#endif //if CUDA_TEST_EXP

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
