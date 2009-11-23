/*
Try accelating splotch with CUDA. July 2009. 
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
#include <math.h>

// includes, project
#include "cuda.h"
#include <cutil_inline.h>

// includes, kernels
#include <splotch_kernel.cu> 
#include "vtimer.h"
#include "splotch_cuda.h"
#include "CuPolicy.h"

///////////////////////////////////////////////////////////////////
//functions defs
//void dump_pr(cu_param_range *pr);
template<typename T> T findParamWithoutChange
//    (paramfile *param, const std::string &key, const T &deflt);
    (paramfile *param,  std::string &key,  T &deflt);
extern "C" void getCuTransformParams(cu_param_transform &p,
    paramfile &params, double campos[3], double lookat[3], double sky[3]);
//usefule functions from splotchutils.h
template<typename T> void get_minmax (T &minv, T &maxv, T val);
////////////////////////////////////////////////////////////////////

//////////////////////////////
//global varibles
/*
float               *d_tmp=0;   //used for debug
float               *d_expTable =0;
CuPolicy            *policy=0;
cu_particle_sim     *d_pd=0;    //device_particle_data
cu_colormap_info    d_colormap_info;    //it contains device pointers
cu_particle_splotch *d_ps_colorize =0; 
cu_exptable_info    d_exp_info; //it contains device pointers
cu_particle_splotch *d_ps_render =0; 
cu_fragment_AeqE    *d_fbuf =0;
cu_color            *d_pic=0;
*/
//////////////////////////////

//////////////////////////////////////////////
//MACROs
#define CLEAR_MEM(p) if(p) {cutilSafeCall(cudaFree(p)); p=0;}
///////////////////////////////////////////////

extern "C" 
void    cu_init(paramfile &params, int devID, cu_gpu_vars* pgv)
{
    //initilize cuda runtime
    //cudaSetDevice( cutGetMaxGflopsDeviceId() );
    cudaSetDevice( devID );
    
    int d;
    cudaGetDevice(&d);
    printf("\nDevice being used %d\n", d);

    unsigned int s;

/*    //d_tmp used for debug
    s =sizeof (float);
    cutilSafeCall( cudaMalloc((void**) &(pgv->d_tmp), s));
    //copy 0.0 to *d_tmp;
    float   f=0.0;
    cutilSafeCall(cudaMemcpy(pgv->d_tmp, &f, sizeof(float),
                              cudaMemcpyHostToDevice) );    
*/

    //Initialize pgv->policy class
    pgv->policy =new CuPolicy(&params);
}

extern "C"
void	cu_range(paramfile &params ,cu_particle_sim* h_pd, 
            unsigned int n, cu_gpu_vars* pgv)
{
   //allocate device memory for particle data
    int s =pgv->policy->GetSizeDPD(n);
#ifdef _DEVICEEMU
    printf("\ndevice_particle_data size:%d\n" ,s);
#endif
    //one more space allocated for the dumb
    cutilSafeCall( cudaMalloc((void**) &pgv->d_pd, s +sizeof(cu_particle_sim)));
    
    //copy particle data to device
    cutilSafeCall(cudaMemcpy(pgv->d_pd, h_pd, s,
                              cudaMemcpyHostToDevice) );    
    //ask for dims from pgv->policy
    dim3    dimGrid, dimBlock;
    pgv->policy->GetDimsRange(n, &dimGrid, &dimBlock);    

    //prepare parameters for stage 1
    cu_param_range  pr;
    int ptypes = params.find<int>("ptypes",1);
    pr.ptypes =ptypes;
    //now collect parameters from configuration
    for(int itype=0;itype<ptypes;itype++)
    {
        pr.log_int[itype] = params.find<bool>("intensity_log"+dataToString(itype),true);
        pr.log_col[itype] = params.find<bool>("color_log"+dataToString(itype),true);
        pr.asinh_col[itype] = params.find<bool>("color_asinh"+dataToString(itype),false);
        pr.col_vector[itype] = params.find<bool>("color_is_vector"+dataToString(itype),false);
        pr.mincol[itype]=1e30;
        pr.maxcol[itype]=-1e30;
        pr.minint[itype]=1e30;
        pr.maxint[itype]=-1e30;
    }
    //allocate memory on device and dump parameters to it
    cu_param_range  *d_pr=0;
    s =sizeof(cu_param_range);
    cutilSafeCall( cudaMalloc((void**) &d_pr, s) );    
    cutilSafeCall(cudaMemcpy(d_pr, &pr, s,
                              cudaMemcpyHostToDevice) );    

    // call device for stage 1
    k_range1<<<dimGrid,dimBlock>>>(d_pr, pgv->d_pd, n);
    
    // copy out particles for min_maxes
    s =pgv->policy->GetSizeDPD(n);
    cutilSafeCall(cudaMemcpy( h_pd, pgv->d_pd,  s,
                              cudaMemcpyDeviceToHost) );    



    //find the min-maxes
    for (int m=0; m<n; m++)
    {
        get_minmax(pr.minint[h_pd[m].type], pr.maxint[h_pd[m].type], h_pd[m].I);
        get_minmax(pr.mincol[h_pd[m].type], pr.maxcol[h_pd[m].type], h_pd[m].C1);
        if (pr.col_vector[h_pd[m].type])
        {
            get_minmax(pr.mincol[h_pd[m].type], pr.maxcol[h_pd[m].type], h_pd[m].C2);
            get_minmax(pr.mincol[h_pd[m].type], pr.maxcol[h_pd[m].type], h_pd[m].C3);
        }
    }
    //for debug: dump_pr(&pr);


    //??copy pr back into device??
    //maybe not needed!

    // call device for stage 2 ptypes times
    // prepare parameters1 first
    float minval_int, maxval_int, minval_col, maxval_col;
    
    for(int itype=0;itype<ptypes;itype++)
    {
        minval_int =findParamWithoutChange<float>(&params,  //in mid of developing only
            "intensity_min"+dataToString(itype),pr.minint[itype]);
        maxval_int = findParamWithoutChange<float>(&params, 
            "intensity_max"+dataToString(itype),pr.maxint[itype]);
        minval_col = findParamWithoutChange<float>(&params, 
            "color_min"+dataToString(itype),pr.mincol[itype]);
        maxval_col = findParamWithoutChange<float>(&params, 
            "color_max"+dataToString(itype),pr.maxcol[itype]);
    //for debug: printf("\n%f, %f, %f, %f\n", minval_int, maxval_int, minval_col, maxval_col);

        k_range2<<<dimGrid, dimBlock>>>(d_pr, pgv->d_pd, n, itype,
            minval_int,maxval_int,minval_col,maxval_col);
    }

    //copy result out to host
    //in mid of development only!!!
    cutilSafeCall(cudaMemcpy(h_pd, pgv->d_pd, s,
                              cudaMemcpyDeviceToHost) );    
    
    //free parameters on device
     CLEAR_MEM((d_pr));

    //pgv->d_pd will be freed in cu_end
}//cu range over

extern "C" void	cu_transform
(paramfile &fparams, unsigned int n,
 double c[3], double l[3], double s[3],
cu_particle_sim* h_pd, cu_gpu_vars* pgv)
{
    //retrieve pamaters for transformaiont first
    cu_param_transform tparams;
    getCuTransformParams(tparams,fparams,c,l,s);

    //arrange memory for the parameters and copy to device
    cu_param_transform  *d_pt;
    int size =sizeof(cu_param_transform);
    cutilSafeCall( cudaMalloc((void**) &d_pt, size) );    
    cutilSafeCall(cudaMemcpy(d_pt, &tparams, size,
                              cudaMemcpyHostToDevice) );    

    //Get block dim and grid dim from pgv->policy object
    dim3    dimGrid, dimBlock;
    pgv->policy->GetDimsRange(n, &dimGrid, &dimBlock);    

    //call device transformation
    k_transform<<<dimGrid,dimBlock>>>(pgv->d_pd, n, d_pt);    

    //free parameters' device memory
     CLEAR_MEM((d_pt));

   //copy result out to host
    //in mid of development only!!!
    size =pgv->policy->GetSizeDPD(n);
    cutilSafeCall(cudaMemcpy(h_pd, pgv->d_pd, size,
                              cudaMemcpyDeviceToHost) );    
 
}

extern "C" void		cu_init_colormap(cu_colormap_info h_info, cu_gpu_vars* pgv)
{
    //allocate memories for colormap and ptype_pionts
    //and dump host data into it
    int size =sizeof(cu_color_map_entry) *h_info.mapSize;    
    cutilSafeCall( cudaMalloc((void**) &pgv->d_colormap_info.map, size));
    cutilSafeCall(cudaMemcpy(pgv->d_colormap_info.map, h_info.map, 
        size, cudaMemcpyHostToDevice) );    
    //type
    size =sizeof(int) *h_info.ptypes;
    cutilSafeCall( cudaMalloc((void**) &pgv->d_colormap_info.ptype_points, size));
    cutilSafeCall(cudaMemcpy(pgv->d_colormap_info.ptype_points, h_info.ptype_points, 
        size, cudaMemcpyHostToDevice) );    
       
    //set fields of global varible pgv->d_colormap_info
    pgv->d_colormap_info.mapSize =h_info.mapSize;
    pgv->d_colormap_info.ptypes  =h_info.ptypes;
   
}

extern "C"	void	cu_colorize(paramfile &params, cu_particle_splotch *h_ps, 
    int n, cu_gpu_vars* pgv)
{
    //fetch parameters for device calling first
    cu_param_colorize   pcolorize;
    pcolorize.res       = params.find<int>("resolution",200);
    pcolorize.ycut0     = params.find<int>("ycut0",0);
    pcolorize.ycut1     = params.find<int>("ycut1",pcolorize.res);
    pcolorize.zmaxval   = params.find<float>("zmax",1.e23);
    pcolorize.zminval   = params.find<float>("zmin",0.0);
    pcolorize.ptypes    = params.find<int>("ptypes",1);

    for(int itype=0; itype<pcolorize.ptypes; itype++)
    {
      pcolorize.brightness[itype] = params.find<double>("brightness"+dataToString(itype),1.);
      pcolorize.grayabsorb[itype] = params.find<float>("gray_absorption"+dataToString(itype),0.2);
      pcolorize.col_vector[itype] = params.find<bool>("color_is_vector"+dataToString(itype),false);
    }
    pcolorize.rfac=1.5;

    //prepare memory for parameters and dump to device
    cu_param_colorize   *d_param_colorize;
    cutilSafeCall( cudaMalloc((void**) &d_param_colorize, sizeof(cu_param_colorize)));
    cutilSafeCall( cudaMemcpy(d_param_colorize, &pcolorize, sizeof(cu_param_colorize), 
        cudaMemcpyHostToDevice));    

    //now prepare memory for d_particle_splotch.
    //one more for dums
    int size =n* sizeof(cu_particle_splotch);
    cutilSafeCall( cudaMalloc((void**) &pgv->d_ps_colorize, size+sizeof(cu_particle_splotch)));

    //fetch grid dim and block dim and call device
    dim3    dimGrid, dimBlock;
    pgv->policy->GetDimsColorize(n, &dimGrid, &dimBlock);    
    k_colorize<<<dimGrid,dimBlock>>>(d_param_colorize, pgv->d_pd, n, pgv->d_ps_colorize,pgv->d_colormap_info);

    //copy the result out
    cutilSafeCall(cudaMemcpy(h_ps, pgv->d_ps_colorize, size, cudaMemcpyDeviceToHost) );    

    //free params memory
    CLEAR_MEM((d_param_colorize));

    //device particle_sim memory can be freed now!

    //particle_splotch memory on device will be freed in cu_end
}

extern "C"	int		cu_get_max_region(cu_gpu_vars* pgv)
{
    if (!pgv->policy) return -1;

    return pgv->policy->GetMaxRegion();
}

extern "C" int cu_get_fbuf_size(cu_gpu_vars* pgv)
{
    if (!pgv->policy) return -1;

    return pgv->policy->GetFBufSize();
}

extern "C"	void	cu_init_exptab(double maxexp, cu_gpu_vars* pgv)
{
    //set common fileds of pgv->d_exp_info
    pgv->d_exp_info.expfac =pgv->d_exp_info.dim2 / maxexp;
    //now make up tab1 and tab2 in host
    float *h_tab1, *h_tab2;
    int dim1 =pgv->d_exp_info.dim1, dim2 =pgv->d_exp_info.dim2;
    h_tab1 =new float[dim1];
    h_tab2 =new float[dim2];
    for (int m=0; m<dim1; ++m)
    {
        h_tab1[m]=exp(m*dim1/pgv->d_exp_info.expfac);
        h_tab2[m]=exp(m/pgv->d_exp_info.expfac);
#ifdef _DEVICEEMU
        if (m==100)
            printf("\ncu_init_exptab %dth: (%f, %f)\n", 
                m, h_tab1[m], h_tab2[m]);
#endif
        
    }

    //allocate device memory and dump
    int size =sizeof(float) *dim1;
    cutilSafeCall( cudaMalloc((void**) &pgv->d_exp_info.tab1, size));
    cutilSafeCall( cudaMemcpy(pgv->d_exp_info.tab1, h_tab1, size, 
        cudaMemcpyHostToDevice));    
    size =sizeof(float) *dim2;
    cutilSafeCall( cudaMalloc((void**) &pgv->d_exp_info.tab2, size));
    cutilSafeCall( cudaMemcpy(pgv->d_exp_info.tab2, h_tab2, size, 
        cudaMemcpyHostToDevice));    
    
    //delete tab1 and tab2 in host
    delete []h_tab1;
    delete []h_tab2;
}

extern "C"	void	cu_prepare_render(cu_particle_splotch *p, 
    int n, cu_gpu_vars* pgv)
{
    //to free some memory that will not be needed
    printf("\nPrepare device rendering...\n");
    //init exp table
    cu_init_exptab(MAX_EXP, pgv);

    //now we just use the space for colorize for particles
//    pgv->d_ps_render =pgv->d_ps_colorize;
    //allocate new memory as it may grow longer after spliting
    CLEAR_MEM(pgv->d_ps_colorize);
    int size = (n+1) *sizeof(cu_particle_splotch);
    cutilSafeCall( cudaMalloc((void**) &pgv->d_ps_render, size));    

    //copy filtered particles into device
    size = n *sizeof(cu_particle_splotch);
    cutilSafeCall(cudaMemcpy(pgv->d_ps_render, p,size, 
        cudaMemcpyHostToDevice) );

    //allocate fragment buffer memory on device
    size =cu_get_fbuf_size( pgv);
    cutilSafeCall( cudaMalloc((void**) &pgv->d_fbuf, size));    
}

extern "C"	void	cu_render1
(int startP, int endP, bool a_eq_e, double grayabsorb, cu_gpu_vars* pgv)
{
    //endP actually exceed the last one to render
    //get dims from pgv->policy object first
    dim3 dimGrid, dimBlock;
    pgv->policy->GetDimsRender(endP-startP, &dimGrid, &dimBlock);

    //call device
    k_render1<<<dimGrid, dimBlock>>>(pgv->d_ps_render, startP, endP, 
        pgv->d_fbuf, a_eq_e, grayabsorb,pgv->d_exp_info);
    
}

extern "C"	void	cu_get_fbuf
(void *h_fbuf, bool a_eq_e, unsigned long n, cu_gpu_vars* pgv)
{
    int size;
    if (a_eq_e)
        size =n* sizeof(cu_fragment_AeqE);
    else
        size =n* sizeof(cu_fragment_AneqE);

    cutilSafeCall(cudaMemcpy(h_fbuf, pgv->d_fbuf,size, 
        cudaMemcpyDeviceToHost) );
}

extern "C"
void	cu_end(cu_gpu_vars* pgv)
{
    // clean up memory
//    CLEAR_MEM((d_tmp));
    CLEAR_MEM((pgv->d_pd));
    CLEAR_MEM((pgv->d_ps_colorize));
    CLEAR_MEM((pgv->d_colormap_info.map));
    CLEAR_MEM((pgv->d_colormap_info.ptype_points));
    CLEAR_MEM((pgv->d_exp_info.tab1));
    CLEAR_MEM((pgv->d_exp_info.tab2));
    CLEAR_MEM((pgv->d_fbuf));
    CLEAR_MEM((pgv->d_pic));

    cudaThreadExit();

    //clear pgv->policy object
    if (pgv->policy)
        delete pgv->policy;
}



template<typename T> T findParamWithoutChange
//    (paramfile *param, const std::string &key, const T &deflt)
    (paramfile *param, std::string &key, T &deflt)
{
    T   value;
    if (param->param_present(key))
    {
        param->findParam(key, value);
        return value;
    }
    else
        return deflt;
}


///////////////////////////////////////////////////////////////////////
////////// NOT USED ////////////////////
#ifdef NOT_USED
extern "C" void cu_copy_particle_sim_to_device
(cu_particle_sim *h_p, int n)
{
    if (!d_pd) return;
    cutilSafeCall(cudaMemcpy(d_pd, h_p, n* sizeof(cu_particle_sim), 
        cudaMemcpyHostToDevice) );  
}

extern "C"  void	cu_post_process(int xres, int yres)
{
    if (!d_pic || !policy) return;

    dim3 dimGrid, dimBlock;
    policy->GetDimsPostProcess(xres, yres, &dimGrid, &dimBlock);
    k_post_process<<<dimGrid, dimBlock>>>(d_pic, xres*yres, d_exp_info);        
}

extern "C"  void	cu_get_pic(cu_color *h_pic, int xres, int yres)
{
    if (!d_pic) return;

    int size =xres* yres * sizeof(cu_color);
    cutilSafeCall(cudaMemcpy(h_pic, d_pic,size, 
        cudaMemcpyDeviceToHost) );  
}

extern "C"  void	cu_init_pic(unsigned int xres, unsigned int yres)
{
    int size = (xres * yres +1)* sizeof(cu_color) ;
    cutilSafeCall( cudaMalloc((void**) &d_pic, size));
}

extern "C"	void	cu_combine(cu_param_combine info)
{
    info.fbuf =d_fbuf;
    info.p  =d_ps_render;
    dim3    dimGrid, dimBlock;
    policy->GetDimsCombine((unsigned int)info.minx, (unsigned int)info.miny, 
        (unsigned int)info.maxx, (unsigned int)info.maxy, &dimGrid, &dimBlock);
    
    //pic must be initilized!
    if (!d_pic) return;

    //call device
    k_combine<<<dimGrid, dimBlock>>>(info.minx, info.miny, info.maxx, info.maxy,
        info.xres, info.yres, d_ps_render, info.pStart, info.pEnd, d_fbuf, d_pic);
}




extern "C"	void	cu_render //one-go render only
(cu_particle_splotch *p, unsigned int n,
 /*int xres, int yres,*/ bool a_eq_e,double grayabsorb)
{
/* moved to cu_prepare_render
    //now we just use the space for colorize for particles
    d_ps_render =d_ps_colorize;
    //allocate new memory as it may grow longer after spliting

    //copy filtered particles into device
    int size = (n+1) *sizeof(cu_particle_splotch);
    cutilSafeCall(cudaMemcpy(d_ps_render, p,size, 
        cudaMemcpyHostToDevice) );

    //allocate fragment buffer memory on device
    size =policy->GetFBufSize(a_eq_e);
    cutilSafeCall( cudaMalloc((void**) &d_fbuf, size));    
*/

    //get dims from policy object first
    dim3 dimGrid, dimBlock;
    policy->GetDimsRender(n, &dimGrid, &dimBlock);

    //call device
    k_render<<<dimGrid, dimBlock>>>(d_ps_render, n, d_fbuf,
        a_eq_e, grayabsorb,d_exp_info);//,xres, yres);

    //fragment buffer will be copied out in cu_get_fbuf();
}

extern "C"	float	cu_get_exp(float arg)
{
    float   *d_result;
    cutilSafeCall( cudaMalloc((void**) &d_result, sizeof(float)));
    
    //in test phase it's a global. after that it's a device func.
    k_get_exp<<<1,1>>>(arg, d_exp_info, d_result);
    
    float   result;
    cutilSafeCall(cudaMemcpy(&result, d_result,sizeof(float), 
        cudaMemcpyDeviceToHost) );    
        
    CLEAR_MEM(d_result);

    return result;
    
}

//cu_get_color is only for debug
extern "C" cu_color	cu_get_color(int ptype, float val)
{
    cu_color   clr, *d_clr;
    cutilSafeCall( cudaMalloc((void**) &d_clr, sizeof(cu_color)));
    
    k_get_color<<<1,1>>>(ptype, val, d_colormap_info, d_clr);
    
    cutilSafeCall(cudaMemcpy(&clr, d_clr,sizeof(cu_color), 
        cudaMemcpyDeviceToHost) );    
        
    CLEAR_MEM((d_clr));

    return clr;
}

void dump_pr(cu_param_range *pr)
{
    printf("\ndump_pr:\n");
    printf("col_vector, log_int,log_col,asinh_col,");
    printf("mincol, maxcol, minint,maxint\n");

    for (int i=0; i<pr->ptypes; i++)
    {
        printf("%d, %d, %d, %d, %f, %f, %f, %f\n",
            pr->col_vector[i], pr->log_int[i], 
            pr->log_col[i],	pr->asinh_col[i], pr->mincol[i],
    		pr->maxcol[i], pr->minint[i], pr->maxint[i]);
    }                     
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
    CLEAR_MEM((d_tmp));
    CLEAR_MEM((d_expTable));
    CLEAR_MEM((d_g_vars));
    CLEAR_MEM((d_p));
    CLEAR_MEM((d_f));
    
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
    CLEAR_MEM((d_result));
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
    CLEAR_MEM((d_expTable));
    CLEAR_MEM((d_tmp));    

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
#endif  //if def NOT_USED
