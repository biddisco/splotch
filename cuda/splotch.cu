/*
Try accelating splotch with CUDA. July 2009. 
*/


// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "cuda.h"
#include <cutil_inline.h>

// includes, kernels
#include <splotch_kernel.cu> 
#include "vtimer.h"
#include "splotch_cuda.h"

//////////////////////////////
//global varibles
float       *d_tmp=0;   //used for debug
float       *d_expTable =0;
//////////////////////////////

extern "C" 
void    cu_init()
{
    //initilize cuda runtime
    cudaSetDevice( cutGetMaxGflopsDeviceId() );

    unsigned int s;

    //d_tmp used for debug
    s =sizeof (float);
    cutilSafeCall( cudaMalloc((void**) &d_tmp, s));
    //copy 0.0 to *d_tmp;
    float   f=0.0;
    cutilSafeCall(cudaMemcpy(d_tmp, &f, sizeof(float),
                              cudaMemcpyHostToDevice) );    
}

extern "C"
void	cu_end()
{
    // clean up memory
    cutilSafeCall(cudaFree(d_tmp));
    cudaThreadExit();
}

/*
extern "C" 
void    cu_initExp(int nExp, float *h_expTable)
{
    //exp table 
    int    s= sizeof(float) *nExp;
    cutilSafeCall( cudaMalloc((void**) &d_expTable, s));        
    cutilSafeCall(cudaMemcpy(d_expTable, h_expTable, s,
                              cudaMemcpyHostToDevice) );    

}

extern "C" 
void    cu_initGVars(G_VARS *h_vars)
{
    //g_vars
    int    s =sizeof(G_VARS);
    cutilSafeCall( cudaMalloc((void**) &d_g_vars, s));
    cutilSafeCall(cudaMemcpy(d_g_vars, h_vars, s,
                              cudaMemcpyHostToDevice) );    

}

extern "C" 
void    cu_initPArray(PARTICLE *h_p)
{
    //particle array
    int    s =sizeof(PARTICLE) *nParticle;
    cutilSafeCall( cudaMalloc((void**) &d_p, s));
    cutilSafeCall(cudaMemcpy(d_p, h_p, s,
                              cudaMemcpyHostToDevice) );    

}

extern "C" 
void    cu_copyPArray2host(PARTICLE *h_p)
{
    //particle array
    int    s =sizeof(PARTICLE) *nParticle;
    cutilSafeCall(cudaMemcpy( h_p, d_p, s,
                              cudaMemcpyDeviceToHost) );    
}

extern "C" 
void    cu_preCalc()
{
    k_preCalc<<<512,512>>>( d_g_vars, d_p);
}

extern "C" 
void    cu_init_frameBuf()
{
    int    s =sizeof(FRAGMENT) *800*800 *64;//is it ok for 64 frame buffers? Yes.
    cutilSafeCall( cudaMalloc((void**) &d_f, s));
    
    dim3 dimGrid(800,800);
    k_initFBuf<<<dimGrid,64>>>(d_f);
}

extern "C"
void    cu_copyFrameBuf2host(FRAGMENT   *h_f)
{
    int    s =sizeof(FRAGMENT) *800*800;
    cutilSafeCall(cudaMemcpy( h_f, d_f, s,
                              cudaMemcpyDeviceToHost) );            
//    int error=cudaMemcpy( h_f, d_f+63*800*800, s,//......
//                cudaMemcpyDeviceToHost);
//    printf("\nerrorCode=%d",error);
}

extern "C"
void    cu_shadeA(int nStart, int nEnd) //do shading due to Plan A
{
    k_shadeA<<<64,1>>> (nStart, nEnd, d_p, d_f, d_g_vars,d_expTable);
}

extern "C"
void	cu_end()
{
    // clean up memory
    cutilSafeCall(cudaFree(d_tmp));
    cutilSafeCall(cudaFree(d_expTable));
    cutilSafeCall(cudaFree(d_g_vars));
    cutilSafeCall(cudaFree(d_p));
    cutilSafeCall(cudaFree(d_f));
    
    cudaThreadExit();
}

extern "C"
void	cu_combineA()
{
    dim3    dimGrid(800,800);//dimGrid(40,40); with 400
    k_combineA<<<dimGrid, 1>>>(d_f);
}
*/

#ifdef CU_DO_TESTS
/*
extern "C"
void	cu_check1()
{
    //check if d_g_var is correctly assigned
    G_VARS  gv;
    cutilSafeCall(cudaMemcpy(&gv, d_g_vars, sizeof(gv),
                              cudaMemcpyDeviceToHost) );
    printf("\ncu_check1: %f,%f,%f,%f,%f,%f,%d,%d,%d",
        gv.rfac, gv.bfak, gv.i00, gv.sigma0, gv.brightness,
        gv.grayabsorb, gv.res, gv.ycut0, gv.ycut1);

    //check if d_p is correctly assigned for computing
    PARTICLE    *p =new PARTICLE[nParticle];
    cutilSafeCall(cudaMemcpy(p, d_p, sizeof(PARTICLE)*nParticle,
                              cudaMemcpyDeviceToHost) );
//    not allowed! cutilSafeCall(cudaMemcpy(&p2, d_p+ (nParticle-1)*sizeof(p1), sizeof(p1),
  //                            cudaMemcpyDeviceToHost) );
//    cutilSafeCall(cudaMemcpy(&p3, d_p+ (nParticle-1)/2*sizeof(p1),sizeof(p1),
//                              cudaMemcpyDeviceToHost) );
    printf("\ncu_check1: %f,%f,%f,%f,%f,%f,%f,%d",
        p[0].x, p[0].y, p[0].z, p[0].r, p[0].ro, p[0].I, p[0].T, p[0].type);
    printf("\ncu_check1: %f,%f,%f,%f,%f,%f,%f,%d",
        p[8888].x, p[8888].y, p[8888].z, p[8888].r, p[8888].ro, p[8888].I, p[8888].T, p[8888].type);
    printf("\ncu_check1: %f,%f,%f,%f,%f,%f,%f,%d",
        p[nParticle-1].x, p[nParticle-1].y, p[nParticle-1].z, p[nParticle-1].r, p[nParticle-1].ro, p[nParticle-1].I, p[nParticle-1].T, p[nParticle-1].type);
    delete p;
    
}


extern "C"
float   cu_testExp(float arg)
{
    float   result;
    //call kernel1 to do xexp( arg);
    k_xexp<<<1,1>>> (10000,-20., arg, d_tmp, d_expTable);
    //copy the result on device back to host
    unsigned int    s =sizeof (float);
//    cutilSafeCall(cudaMemcpy(&result, d_tmp, s,
//                              cudaMemcpyDeviceToHost) );    //highly time consuming
    
    //return    it to caller
    return result;
}

extern "C"
double cu_TestDouble()
{
    float  *d_result, result;
    int s =sizeof(result);
    
    cutilSafeCall( cudaMalloc((void**) &d_result, s));
    k_testDouble<<<1,1>>> (d_result);
    cutilSafeCall(cudaMemcpy(&result, d_result, s,
                              cudaMemcpyDeviceToHost) );
    cutilSafeCall(cudaFree(d_result));
    return result;    
}



extern "C"
float cu_test1(int n, float  *table, float arg)
{
    float result =1;

    //allocate space on device for exp table 
    unsigned int s =sizeof(float) *n;
    float   *d_expTable =0;
    cutilSafeCall( cudaMalloc((void**) &d_expTable, s));        

    //copy exp table from host to device
    cutilSafeCall(cudaMemcpy(d_expTable, table, s,
                              cudaMemcpyHostToDevice) );    
    

    //verify things
    //allcate a host mem for copying device data out
    float *h_temp =(float*)malloc(s);
    //copy device exp table to this temp memory, from device to host
    cutilSafeCall(cudaMemcpy(h_temp, d_expTable, s,
                              cudaMemcpyDeviceToHost) );    
    //compare the two tables
    for (int i=0; i<n; i++)
    {
        if ( h_temp[i] != table[i] )
        {
            result =-1;
            break;
        }
        if ( i==0 || i==n/2 || i==n-1) //pick some to check
            printf("\n %f == %f", h_temp[i], table[i]);
    }

    //test calling a kernel function to retrieve a table element
    //allocate the memory for the single float result on device
    float   *d_tmp;
    s =sizeof (float);
    cutilSafeCall( cudaMalloc((void**) &d_tmp, s));

    //call kernel to get d_expTable[10] to d_tmp
    test_kernel <<<1,1>>>( d_expTable, 10, d_tmp);       
    //copy the result on device back to host
    float   a;
    cutilSafeCall(cudaMemcpy(&a, d_tmp, s,
                              cudaMemcpyDeviceToHost) );    
    //and compare
    if ( table[10] != a)
        result =-1;

    //call kernel1 to do xexp( arg);
    k_xexp<<<1,1>>> (10000,-20., arg, d_tmp, d_expTable);
    //copy the result on device back to host
    cutilSafeCall(cudaMemcpy(&a, d_tmp, s,
                              cudaMemcpyDeviceToHost) );    
    //return    it to caller
    result =a;

    //free the memory on host
    free(h_temp);

    //free the device memory
    cutilSafeCall(cudaFree(d_expTable));
    cutilSafeCall(cudaFree(d_tmp));    

    return result;
}


extern "C"
void cu_test2(int *l)
{
    VTimer t;

    int    *d_l;
    int sz =sizeof(int) *65535*512;
    

    cutilSafeCall( cudaMalloc((void**) &d_l, sz));       
    
    t.start();
    test_kernel2<<<65535,512>>>(d_l);
	t.stop();
	printf("\ntime1=%f", t.getTime());
    
    int err= cudaMemcpy(l, d_l, sz, cudaMemcpyDeviceToHost) ;    

    cudaFree(d_l);  
}
*/
#endif //CU_DO_TESTS