#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
#define MAX_P_TYPE 8//('XXXX','TEMP','U','RHO','MACH','DTEG','DISS','VEL')
                                        //in mid of developing only
#define MAX_EXP -20.0
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable

typedef struct 
  {
  float r,g,b;
  }cu_color;

 typedef struct 
  {
  cu_color e;
  float x,y,z,r,I;
  unsigned short type;
  bool active;
  
  } cu_particle_sim;

typedef struct 
  {
	cu_color e;
  float x,y,z,r,I;
  unsigned short type;
  bool active;
  
  } particle_sim2;



typedef struct 
  {
  float x,y,r,I;
  int type;
  cu_color e;
  bool isValid;
  unsigned short minx, miny, maxx, maxy;
  unsigned long posInFragBuf;
  } cu_particle_splotch;

typedef struct 
  {
  float p[12];
  bool  projection;
  int   xres, yres;
  float fovfct, dist, xfac;
  float minrad_pix;
  int ptypes;
  float zmaxval, zminval;
  bool col_vector[MAX_P_TYPE];
  float brightness[MAX_P_TYPE];
  float rfac;
  }cu_param;




typedef struct 
  {
  float val;
  cu_color color;
  }cu_color_map_entry;

typedef struct 
  {
  cu_color_map_entry *map;
  int mapSize;
  int *ptype_points;
  int ptypes;
  }cu_colormap_info;

typedef struct 
  {
  float aR, aG, aB;
  } cu_fragment_AeqE;

typedef struct 
  {
  float aR, aG, aB;
  float qR, qG, qB;
  } cu_fragment_AneqE;


//MACROs
#define Pi 3.14159265358979323846264338327950288

#define MAXSIZE 1000



//fetch a color from color table on device
 cu_color get_color(int ptype, float val, int mapSize, int ptypes, cu_color_map_entry *dmap_in, int *ptype_points_in) 
  {
  	//__local int map_size; //not working on AMD
	 int map_size;
	//__local int map_ptypes;not working on AMD
	 int map_ptypes;
int ptype_points[10];
  cu_color_map_entry dmap[MAXSIZE];
for (int i=0;i<ptypes;i++) ptype_points[i]=ptype_points_in[i];
for (int i=0;i<mapSize;i++) dmap[i]=dmap_in[i];

  map_size = mapSize;
  map_ptypes = ptypes;
  //first find the right entry for this ptype
  int     start, end;

  start =ptype_points[ptype];
  if ( ptype == map_ptypes-1)//the last type
    end =map_size-1;
  else
    end =ptype_points[ptype+1]-1;

  //search the section of this type to find the val
  int i=start;
  while ((val>dmap[i+1].val) && (i<end)) ++i;

  const float fract = (val-dmap[i].val)/(dmap[i+1].val-dmap[i].val);
  cu_color clr1=dmap[i].color, clr2=dmap[i+1].color;
  cu_color        clr;
  clr.r =clr1.r + fract*(clr2.r-clr1.r);
  clr.g =clr1.g + fract*(clr2.g-clr1.g);
  clr.b =clr1.b + fract*(clr2.b-clr1.b);

  return clr;
  }


//device render function k_render1
__kernel void k_render1
  (__global cu_particle_splotch *p, int nP,__global 
  cu_fragment_AeqE *buf,  float grayabsorb, int mapSize, int types,__global cu_param *dparams1,__global cu_color_map_entry *dmap_in,__global int *ptype_points_in)
  {
  
  int ptype_points[10];
  cu_color_map_entry dmap[MAXSIZE];
for (int i=0;i<types;i++) ptype_points[i]=ptype_points_in[i];
for (int i=0;i<mapSize;i++) dmap[i]=dmap_in[i];
cu_param dparams = dparams1[0];
  //first get the index m of this thread
  int m;
//__constant__ cu_param dparams;
  //m =get_group_id(0)*get_local_size(0)+get_local_id(0);
  m = get_global_id(0);
  if (m<nP) {

  // coloring
	int ptype = p[m].type;
	float col1=p[m].e.r,col2=p[m].e.g,col3=p[m].e.b;
	clamp (0.0000001,0.9999999,col1);
	if (dparams1->col_vector[ptype])
	{
		col2 = clamp (0.0000001,0.9999999,col2);
		col3 = clamp (0.0000001,0.9999999,col3);
	}
	float intensity=p[m].I;
	intensity = clamp (0.0000001,0.9999999,intensity);
	intensity *= dparams1->brightness[ptype];

	cu_color e;
	if (dparams1->col_vector[ptype]) // color from file
	{
		e.r=col1*intensity;
		e.g=col2*intensity;
		e.b=col3*intensity;
	}
	else // get color, associated from physical quantity contained in e.r, from lookup table
	{
		//first find the right entry for this ptype
		if (ptype<types)
		{
			e = get_color(ptype, col1, mapSize, types,dmap,ptype_points);
			e.r *= intensity;
			e.g *= intensity;
			e.b *= intensity;
		}
		else
		{	e.r =e.g =e.b =0.0;}
	}

	//make fbuf the right type

	//now do the rendering
	const float powtmp = pow(Pi,1./3.);
	const float sigma0 = powtmp/sqrt(2*Pi);

	const float r = p[m].r;
	const float radsq = 2.25*r*r;
	const float stp = -0.5/(r*r*sigma0*sigma0);

	cu_color q;//e=p[m].e;
	
	const float intens = -0.5/(2*sqrt(Pi)*powtmp);
	e.r*=intens; e.g*=intens; e.b*=intens;

	const float posx=p[m].x, posy=p[m].y;
	unsigned int fpos =p[m].posInFragBuf;

	
		for (int x=p[m].minx; x<p[m].maxx; ++x)
		{
			float dxsq=(x-posx)*(x-posx);
			for (int y=p[m].miny; y<p[m].maxy; ++y) 
			{
				float dsq = (y-posy)*(y-posy) + dxsq;
				if (dsq<radsq)
				{
					float att = pow((float)2.71828,(stp*dsq));// = expf(stp*dsq); //$$$$$$$$$$$$$$$$$$$$$$$$$$
				buf[fpos].aR = att*e.r;
				buf[fpos].aG = att*e.g;
				buf[fpos].aB = att*e.b;
				}
				else
				{
			buf[fpos].aR =0.0;
			buf[fpos].aG =0.0;
			buf[fpos].aB =0.0;
				}
				//for each (x,y)
				fpos++;
			} //y
		} //x
	
	
}

}

__kernel void k_render2
  (__global cu_particle_splotch *p, int nP,__global 
  cu_fragment_AneqE *buf,  float grayabsorb, int mapSize, int types,__global cu_param *dparams1,__global cu_color_map_entry *dmap_in,__global int *ptype_points_in)
  {
  
  int ptype_points[10];
  cu_color_map_entry dmap[MAXSIZE];
for (int i=0;i<types;i++) ptype_points[i]=ptype_points_in[i];
for (int i=0;i<mapSize;i++) dmap[i]=dmap_in[i];
cu_param dparams = dparams1[0];
  //first get the index m of this thread
  int m;
//__constant__ cu_param dparams;
 // m =get_group_id(0)*get_local_size(0)+get_local_id(0);
  m = get_global_id(0);
  if (m<nP) {

  // coloring
	int ptype = p[m].type;
	float col1=p[m].e.r,col2=p[m].e.g,col3=p[m].e.b;
	clamp (0.0000001,0.9999999,col1);
	if (dparams1->col_vector[ptype])
	{
		col2 = clamp (0.0000001,0.9999999,col2);
		col3 = clamp (0.0000001,0.9999999,col3);
	}
	float intensity=p[m].I;
	intensity = clamp (0.0000001,0.9999999,intensity);
	intensity *= dparams1->brightness[ptype];

	cu_color e;
	if (dparams1->col_vector[ptype]) // color from file
	{
		e.r=col1*intensity;
		e.g=col2*intensity;
		e.b=col3*intensity;
	}
	else // get color, associated from physical quantity contained in e.r, from lookup table
	{
		//first find the right entry for this ptype
		if (ptype<types)
		{
			e = get_color(ptype, col1, mapSize, types,dmap,ptype_points);
			e.r *= intensity;
			e.g *= intensity;
			e.b *= intensity;
		}
		else
		{	e.r =e.g =e.b =0.0;}
	}

	//make fbuf the right type
	

	//now do the rendering
	const float powtmp = pow(Pi,1./3.);
	const float sigma0 = powtmp/sqrt(2*Pi);

	const float r = p[m].r;
	const float radsq = 2.25*r*r;
	const float stp = -0.5/(r*r*sigma0*sigma0);

	cu_color q;//e=p[m].e;
	
		q.r = e.r/(e.r+grayabsorb);
		q.g = e.g/(e.g+grayabsorb);
		q.b = e.b/(e.b+grayabsorb);
	
	const float intens = -0.5/(2*sqrt(Pi)*powtmp);
	e.r*=intens; e.g*=intens; e.b*=intens;

	const float posx=p[m].x, posy=p[m].y;
	unsigned int fpos =p[m].posInFragBuf;

	
	
		for (int x=p[m].minx; x<p[m].maxx; ++x)
		{
			float dxsq=(x-posx)*(x-posx);
			for (int y=p[m].miny; y<p[m].maxy; ++y)
			{
				float dsq = (y-posy)*(y-posy) + dxsq;
				if (dsq<radsq)
				{
					float att =pow((float)2.71828,(stp*dsq));////////////= __expf(stp*dsq);
					float expm1 = pow((float)2.71828,(att*e.r))-1.0;
					//expm1 ;//////////////////=__expf(att*e.r)-1.0;
			buf[fpos].aR = expm1;
			buf[fpos].qR = q.r;
					expm1= pow((float)2.71828,(att*e.g))-1.0;/////////////// =__expf(att*e.g)-1.0;
			buf[fpos].aG = expm1;
			buf[fpos].qG = q.g;
					expm1= pow((float)2.71828,(att*e.b))-1.0;//////////// =__expf(att*e.b)-1.0;
			buf[fpos].aB = expm1;
			buf[fpos].qB = q.b;
				}
				else
				{
			buf[fpos].aR =0.0;
			buf[fpos].aG =0.0;
			buf[fpos].aB =0.0;
			buf[fpos].qR =1.0;
			buf[fpos].qG =1.0;
			buf[fpos].qB =1.0;
				}
				//for each (x,y)
				fpos++;
			} //y
		} //x
	}


}

//Transform by kernel
__kernel void k_transform
  (__global particle_sim2 *p,__global  cu_particle_splotch *p2,int n, __global cu_param *dparams1)
  {
  

  //first get the index m of this thread
  
  //int m=get_group_id(0)*get_local_size(0)+get_local_id(0);
  int m=get_global_id(0);
  if (m <=n) {

  //now do x,y,z
  float x,y,z;
  x =p[m].x*dparams1->p[0] + p[m].y*dparams1->p[1] + p[m].z*dparams1->p[2] + dparams1->p[3];
  y =p[m].x*dparams1->p[4] + p[m].y*dparams1->p[5] + p[m].z*dparams1->p[6] + dparams1->p[7];
  z =p[m].x*dparams1->p[8] + p[m].y*dparams1->p[9] + p[m].z*dparams1->p[10]+ dparams1->p[11];

  //do r
  float xfac = dparams1->xfac;
  const float   res2 = 0.5*dparams1->xres;
  const float   ycorr = .5f*(dparams1->yres-dparams1->xres);
  if (!dparams1->projection)
    {
    x = res2 * (x+dparams1->fovfct*dparams1->dist)*xfac;
    y = res2 * (y+dparams1->fovfct*dparams1->dist)*xfac + ycorr;
    }
  else
    {
    xfac=1./(dparams1->fovfct*z);
    x = res2 * (x+dparams1->fovfct*z)*xfac;
    y = res2 * (y+dparams1->fovfct*z)*xfac + ycorr;
    }

  float r = p[m].r;
  p[m].I /= r;
  r *= res2*xfac;

  const float rfac= sqrt(r*r + 0.25*dparams1->minrad_pix*dparams1->minrad_pix)/r;
  r *= rfac;
  p2[m].I = p[m].I/rfac;

  p2[m].isValid = false;

  // compute region occupied by the partile
  const float rfacr=dparams1->rfac*r;
  int minx=(int)(x-rfacr+1);
  if (minx>=dparams1->xres) return;
  minx=max(minx,0);

  int maxx=(int)(x+rfacr+1);
  if (maxx<=0) return;
  maxx=min(maxx,dparams1->xres);
  if (minx>=maxx) return;

  int miny=(int)(y-rfacr+1);
  if (miny>=dparams1->yres) return;
  miny=max(miny,0);

  int maxy=(int)(y+rfacr+1);
  if (maxy<=0) return;
  maxy=min(maxy,dparams1->yres);
  if (miny>=maxy) return;

  p2[m].minx =minx;  p2[m].miny =miny;
  p2[m].maxx =maxx;  p2[m].maxy =maxy;

  p2[m].isValid = true;
  p2[m].x = x;
  p2[m].y = y;
  p2[m].r = r;
  p2[m].e.r =  p[m].e.r;
  p2[m].e.g =  p[m].e.g;
  p2[m].e.b =  p[m].e.b;
  p2[m].type = p[m].type;

  }
}

