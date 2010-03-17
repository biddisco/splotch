#include "cuda/splotchutils_cuda.h"

using namespace std;

extern "C"
void getCuTransformParams(cu_param_transform &para_trans,
paramfile &params, double c[3], double l[3], double s[3])
  {
  vec3 campos(c[0], c[1], c[2]),
       lookat(l[0], l[1], l[2]),
       sky   (s[0], s[1], s[2]);

  int res = params.find<int>("resolution",200);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);
  float64 xfac=0.0, dist=0.0;

  sky.Normalize();
  vec3 zaxis = (lookat-campos).Norm();
  vec3 xaxis = crossprod (sky,zaxis).Norm();
  vec3 yaxis = crossprod (zaxis,xaxis);
  TRANSFORM trans;
  trans.Make_General_Transform
        (TRANSMAT(xaxis.x,xaxis.y,xaxis.z,
                  yaxis.x,yaxis.y,yaxis.z,
                  zaxis.x,zaxis.y,zaxis.z,
                  0,0,0));
  trans.Invert();
  TRANSFORM trans2;
  trans2.Make_Translation_Transform(-campos);
  trans2.Add_Transform(trans);
  trans=trans2;
  bool projection = params.find<bool>("projection",true);

  if (!projection)
    {
    dist= (campos-lookat).Length();
    xfac=1./(fovfct*dist);
    cout << " Field of fiew: " << 1./xfac*2. << endl;
    }

  bool minhsmlpixel = params.find<bool>("minhsmlpixel",false);

  //retrieve the parameters for tansformation
  for (int i=0; i<12; i++)
    para_trans.p[i] =trans.Matrix().p[i];
  para_trans.projection=projection;
  para_trans.res=res;
  para_trans.fovfct=fovfct;
  para_trans.dist=dist;
  para_trans.xfac=xfac;
  para_trans.minhsmlpixel=minhsmlpixel;
  }
